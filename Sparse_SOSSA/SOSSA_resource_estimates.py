import numpy as np
import cvxpy as cp
import mosek


class ResourceEstimates:
    """

    Spin-explicit / spin-free SOS SDP (Appendix G of arXiv:2502.15882v1).

    Variables:
        S_G = [[G, Gp],[Gpp, Gppp]] ∈ S_+^{2P}
          with P = M^2 (ordered pairs) or P = M(M+1)/2 (triangular pairs).
        If spin_explicit:
            D_up, D_dn, Q_up, Q_dn ∈ S_+^{M}
        else:
            D, Q ∈ S_+^{M}

    Constraints:
      (G10)  for all i<=k, j<=l  (triangular)   OR   for all i,k,j,l (ordered):
        G[(k,i),(l,j)] - Gp[(k,i),(j,l)] - Gpp[(i,k),(l,j)] + Gppp[(i,k),(j,l)]
        = 0.5 * h2[i,k,j,l]

      (G11)  for each (i,j):
        2 * Σ_p [Gp[(j,i),(p,p)] + Gpp[(p,p),(j,i)] - Gppp[(p,p),(i,j)] - Gppp[(i,j),(p,p)]]
        +  (D_up - Q_up)[j,i] = h1[i,j]
        and similarly for ↓, or spin-free version:
        2 * Σ_p [ ... ] + (D[i,j] - Q[j,i]) = h1[i,j]

      (G12) Objective (minimize):
        ESOS = 4 * Σ_{p,q} Gppp[(p,p),(q,q)] + tr(Q_up) + tr(Q_dn)      (spin-explicit)
        ESOS = 4 * Σ_{p,q} Gppp[(p,p),(q,q)] + 2 * tr(Q)                (spin-free)

    Variables:
        S_G = [[G, Gp], [Gpp, Gppp]] in PSD cone, shape (2P, 2P).

        If spin_explicit:
            D_up, D_dn, Q_up, Q_dn in PSD cone, shape (M, M).

        Else:
            D, Q in PSD cone, shape (M, M).

    The SDP variables are uppercase Gram/PSD matrices. The actual SOS-generator
    coefficients appearing in O_alpha are obtained after PSD factorization:

        S_G = Y_SF^T Y_SF,
        D   = Y_D^T  Y_D,
        Q   = Y_Q^T  Y_Q.

    Each row of these factor matrices corresponds to one SOS generator. This
    class estimates d_sparse by counting nonzero entries of the lowercase
    generator coefficients:

        g, gbar, d, q

    after thresholding. This is the sparse support estimate appropriate for
    the O_alpha-based SOSSA block-encoding model.

    Returned:
        ESOS, lambda_sqrt2, Lambda, coeffs, diagnostics
    """

    def __init__(
        self,
        norb,
        h1,
        h2,
        *,
        verbose=False,
        solver="SCS",
        max_iterations=50000,
        spin_explicit=True,
        enforce_spin_free=False,
        use_triu_pairs=True,
        h2_abs_threshold=0.0,
        support_threshold=1e-8,
    ):
        self.M = int(norb)  # spatial orbitals
        self.h1 = np.asarray(h1, dtype=float)
        self.h2 = np.asarray(h2, dtype=float)

        assert self.h1.shape == (self.M, self.M), "h1 must be shape (M, M)."
        assert self.h2.shape == (self.M, self.M, self.M, self.M), (
            "h2 must be shape (M, M, M, M)."
        )

        self.verbose = bool(verbose)
        self.solver = solver
        self.max_iters = int(max_iterations)
        self.spin_explicit = bool(spin_explicit)
        self.enforce_spin_free = bool(enforce_spin_free)

        self.use_triu_pairs = bool(use_triu_pairs)

        if self.use_triu_pairs:
            self.P = self.M * (self.M + 1) // 2

            def pidx(i, k):
                if k < i:
                    i, k = k, i
                return i * self.M - (i * (i - 1)) // 2 + (k - i)

            self.pidx = pidx

        else:
            self.P = self.M * self.M
            self.pidx = lambda a, b: a * self.M + b

        # Threshold used while building G10 constraints.
        self.h2_thresh = float(h2_abs_threshold)

        # Threshold used after solving to count nonzero entries of lowercase
        # g, gbar, d, q generator coefficients.
        self.support_threshold = float(support_threshold)

        self._last_num_g10_dropped = 0
        self._last_solution = None
        self._resource_configured = False

    def _diag_pair_indices(self):
        return np.array([self.pidx(p, p) for p in range(self.M)], dtype=int)

    def _build_problem(self):
        M, P = self.M, self.P

        S_G = cp.Variable((2 * P, 2 * P), PSD=True)

        def blk(A, r0, c0, n):
            return A[r0 : r0 + n, c0 : c0 + n]

        G = blk(S_G, 0, 0, P)
        Gp = blk(S_G, 0, P, P)
        Gpp = blk(S_G, P, 0, P)
        Gppp = blk(S_G, P, P, P)

        cons = []

        if self.spin_explicit:
            D_up = cp.Variable((M, M), PSD=True)
            D_dn = cp.Variable((M, M), PSD=True)
            Q_up = cp.Variable((M, M), PSD=True)
            Q_dn = cp.Variable((M, M), PSD=True)

            if self.enforce_spin_free:
                cons += [D_up == D_dn, Q_up == Q_dn]

        else:
            D = cp.Variable((M, M), PSD=True)
            Q = cp.Variable((M, M), PSD=True)

        pidx = self.pidx
        PP = self._diag_pair_indices()

        # ------------------------------------------------------------------
        # G10 two-body constraints
        # ------------------------------------------------------------------
        dropped = 0

        if self.use_triu_pairs:
            i_pairs = [(i, k) for i in range(M) for k in range(i, M)]
            j_pairs = [(j, l) for j in range(M) for l in range(j, M)]
        else:
            i_pairs = [(i, k) for i in range(M) for k in range(M)]
            j_pairs = [(j, l) for j in range(M) for l in range(M)]

        for (i, k) in i_pairs:
            for (j, l) in j_pairs:
                val = 0.5 * self.h2[i, k, j, l]

                if self.h2_thresh > 0.0 and abs(val) < self.h2_thresh:
                    dropped += 1
                    continue

                rG, cG = pidx(k, i), pidx(l, j)
                rGp, cGp = pidx(k, i), pidx(j, l)
                rGpp, cGpp = pidx(i, k), pidx(l, j)
                rGppp, cGppp = pidx(i, k), pidx(j, l)

                cons.append(
                    G[rG, cG]
                    - Gp[rGp, cGp]
                    - Gpp[rGpp, cGpp]
                    + Gppp[rGppp, cGppp]
                    == val
                )

        self._last_num_g10_dropped = dropped

        # ------------------------------------------------------------------
        # G11 one-body constraints
        # ------------------------------------------------------------------
        for i in range(M):
            for j in range(M):
                rij = pidx(j, i)
                i_j = pidx(i, j)

                term = (
                    cp.sum(Gp[rij, PP])
                    + cp.sum(Gpp[PP, rij])
                    - cp.sum(Gppp[PP, i_j])
                    - cp.sum(Gppp[i_j, PP])
                )

                if self.spin_explicit:
                    cons += [
                        2.0 * term + (D_up[j, i] - Q_up[j, i]) == self.h1[i, j],
                        2.0 * term + (D_dn[j, i] - Q_dn[j, i]) == self.h1[i, j],
                    ]
                else:
                    cons += [
                        2.0 * term + (D[i, j] - Q[j, i]) == self.h1[i, j]
                    ]

        # ------------------------------------------------------------------
        # G12 objective
        # ------------------------------------------------------------------
        diag_idx = self._diag_pair_indices()
        sum_Gppp_diag = cp.sum(Gppp[diag_idx, :][:, diag_idx])

        if self.spin_explicit:
            ESOS = 4.0 * sum_Gppp_diag + cp.trace(Q_up) + cp.trace(Q_dn)
        else:
            ESOS = 4.0 * sum_Gppp_diag + 2.0 * cp.trace(Q)

        prob = cp.Problem(cp.Minimize(ESOS), cons)

        if self.spin_explicit:
            return prob, S_G, (D_up, D_dn), (Q_up, Q_dn)

        return prob, S_G, D, Q

    def _unflatten_pair(self, v):
        """
        Convert a pair-vector into an M x M matrix of lowercase generator
        coefficients.

        If use_triu_pairs=True, this symmetrizes the triangular storage.
        If use_triu_pairs=False, this reshapes the ordered-pair vector.
        """
        M = self.M

        if self.use_triu_pairs:
            A = np.zeros((M, M))
            idx = 0

            for i in range(M):
                for k in range(i, M):
                    A[i, k] = v[idx]
                    A[k, i] = v[idx]
                    idx += 1

            return A

        return v.reshape(M, M, order="C")

    @staticmethod
    def _psd_factor(A, rtol=1e-7, atol=1e-7):
        """
        Factor a PSD matrix A into rows Y such that approximately:

            A = Y.T @ Y.

        Each row of Y is one SOS-generator coefficient vector.
        """
        A = 0.5 * (A + A.T)
        w, V = np.linalg.eigh(A)
        w = np.clip(w, 0.0, None)

        if w.size == 0:
            return np.zeros((0, A.shape[0])), 0

        wmax = float(np.max(w))
        if wmax <= 0.0:
            return np.zeros((0, A.shape[0])), 0

        keep = w > max(atol, rtol * wmax)

        if not np.any(keep):
            return np.zeros((0, A.shape[0])), 0

        Y = np.sqrt(w[keep])[:, None] * V[:, keep].T
        return Y, int(np.count_nonzero(keep))

    def _estimate_sparse_d_from_generator_coeffs(self, coeffs, threshold=None):
        """
        Estimate d_sparse from lowercase SOS-generator coefficients:

            g, gbar, d, q.

        This is the desired sparse-support estimate for O_alpha operators
        appearing in G1, G2, and G3.

        Spin-free convention:
            coeffs["d_up"] stores the spin-free d factors.
            coeffs["q_up"] stores the spin-free q factors.
            coeffs["d_dn"] and coeffs["q_dn"] are empty placeholders.

        Spin-explicit convention:
            d_up, d_dn, q_up, q_dn are counted separately.

        Formula:
            d_sparse =
                identity_terms
              + nnz(g)
              + nnz(gbar)
              + nnz(d)
              + nnz(q).
        """
        if threshold is None:
            threshold = self.support_threshold

        threshold = float(threshold)

        g = np.asarray(coeffs["g"])
        gbar = np.asarray(coeffs["gbar"])
        d_up = np.asarray(coeffs["d_up"])
        d_dn = np.asarray(coeffs["d_dn"])
        q_up = np.asarray(coeffs["q_up"])
        q_dn = np.asarray(coeffs["q_dn"])

        nnz_g = int(np.count_nonzero(np.abs(g) > threshold))
        nnz_gbar = int(np.count_nonzero(np.abs(gbar) > threshold))

        nnz_d_up = int(np.count_nonzero(np.abs(d_up) > threshold))
        nnz_d_dn = int(np.count_nonzero(np.abs(d_dn) > threshold))
        nnz_q_up = int(np.count_nonzero(np.abs(q_up) > threshold))
        nnz_q_dn = int(np.count_nonzero(np.abs(q_dn) > threshold))

        if self.spin_explicit:
            nnz_d = nnz_d_up + nnz_d_dn
            nnz_q = nnz_q_up + nnz_q_dn
        else:
            nnz_d = nnz_d_up
            nnz_q = nnz_q_up

        # generator coefficients only.

        d_sparse = int(nnz_g + nnz_gbar + nnz_d + nnz_q)

        return {
            "support_model": "generator_coefficients_g_gbar_d_q",
            "threshold": threshold,
            "d_sparse": d_sparse,

            "nnz_g": nnz_g,
            "nnz_gbar": nnz_gbar,
            "nnz_d": nnz_d,
            "nnz_q": nnz_q,

            "nnz_d_up": nnz_d_up,
            "nnz_d_dn": nnz_d_dn,
            "nnz_q_up": nnz_q_up,
            "nnz_q_dn": nnz_q_dn,

            "rank_g_gbar": int(g.shape[0]) if g.ndim >= 1 else 0,
            "rank_d_up": int(d_up.shape[0]) if d_up.ndim >= 1 else 0,
            "rank_d_dn": int(d_dn.shape[0]) if d_dn.ndim >= 1 else 0,
            "rank_q_up": int(q_up.shape[0]) if q_up.ndim >= 1 else 0,
            "rank_q_dn": int(q_dn.shape[0]) if q_dn.ndim >= 1 else 0,

            "g_shape": tuple(g.shape),
            "gbar_shape": tuple(gbar.shape),
            "d_up_shape": tuple(d_up.shape),
            "d_dn_shape": tuple(d_dn.shape),
            "q_up_shape": tuple(q_up.shape),
            "q_dn_shape": tuple(q_dn.shape),
        }

    def solve(self):
        prob, S_G, Dblk, Qblk = self._build_problem()

        kwargs = dict(verbose=self.verbose)
        s = self.solver.upper()

        if s == "MOSEK":
            kwargs.update(
                dict(
                    solver=cp.MOSEK,
                    mosek_params={
                        "MSK_DPAR_INTPNT_CO_TOL_PFEAS": 1e-4,
                        "MSK_DPAR_INTPNT_CO_TOL_DFEAS": 1e-4,
                        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP": 1e-4,
                        "MSK_IPAR_INTPNT_MAX_NUM_COR": -1,
                    },
                )
            )

        elif s == "SCS":
            kwargs.update(
                dict(
                    solver=cp.SCS,
                    max_iters=self.max_iters,
                    eps=1e-4,
                    eps_abs=1e-4,
                    eps_rel=1e-4,
                    alpha=1.8,
                    normalize=True,
                    acceleration_lookback=30,
                    use_indirect=False,
                )
            )

        else:
            kwargs.update(dict(solver=self.solver))

        prob.solve(**kwargs)
        status = prob.status

        M, P = self.M, self.P

        empty = dict(
            g=np.zeros((0, M, M)),
            gbar=np.zeros((0, M, M)),
            d_up=np.zeros((0, M)),
            d_dn=np.zeros((0, M)),
            q_up=np.zeros((0, M)),
            q_dn=np.zeros((0, M)),
        )

        base = dict(
            ESOS=None,
            lambda_sqrt2=None,
            Lambda=None,
            coeffs=empty,
            diagnostics=dict(
                status=status,
                message=None,
                rank_S=0,
                rank_D_up=0,
                rank_D_dn=0,
                rank_Q_up=0,
                rank_Q_dn=0,
                residual_one_e_maxabs=None,
                g10_rows_dropped=int(self._last_num_g10_dropped),
            ),
            total_rank=0,
            sos_support=None,
        )

        if status not in ("optimal", "optimal_inaccurate"):
            base["diagnostics"]["message"] = "Solve failed or infeasible."
            self._last_solution = base
            return base

        # ------------------------------------------------------------------
        # Extract solved SDP blocks.
        # ------------------------------------------------------------------
        Sval = 0.5 * (S_G.value + S_G.value.T)

        G = Sval[0:P, 0:P]
        Gp = Sval[0:P, P : 2 * P]
        Gpp = Sval[P : 2 * P, 0:P]
        Gppp = Sval[P : 2 * P, P : 2 * P]

        if self.spin_explicit:
            D_up, D_dn = Dblk
            Q_up, Q_dn = Qblk

            D_up = 0.5 * (D_up.value + D_up.value.T)
            D_dn = 0.5 * (D_dn.value + D_dn.value.T)
            Q_up = 0.5 * (Q_up.value + Q_up.value.T)
            Q_dn = 0.5 * (Q_dn.value + Q_dn.value.T)

            lambda_sqrt2 = (
                np.trace(G)
                + np.trace(Gppp)
                + np.trace(D_up)
                + np.trace(D_dn)
                + np.trace(Q_up)
                + np.trace(Q_dn)
            )

        else:
            D = 0.5 * (Dblk.value + Dblk.value.T)
            Q = 0.5 * (Qblk.value + Qblk.value.T)

            lambda_sqrt2 = (
                np.trace(G)
                + np.trace(Gppp)
                + 2.0 * np.trace(D)
                + 2.0 * np.trace(Q)
            )

        ESOS = float(prob.value)
        Lambda_eq9 = 0.5 * lambda_sqrt2

        # ------------------------------------------------------------------
        # Factor S_G into lowercase g and gbar coefficients.
        #
        # S_G is the Gram matrix for the spin-free generator:
        #
        #   O_SF = sum_ij g_ij E_ji + gbar_ij hole_rotation_ji
        #
        # Each row of Y_SF gives one alpha generator.
        # ------------------------------------------------------------------
        Y_SF, rank_S = self._psd_factor(Sval, rtol=1e-7, atol=1e-6)

        if rank_S > 0:
            Y_g = Y_SF[:, :P]
            Y_gbar = Y_SF[:, P:]
        else:
            Y_g = np.zeros((0, P))
            Y_gbar = np.zeros((0, P))

        g_list = [self._unflatten_pair(Y_g[a]) for a in range(Y_g.shape[0])]
        gbar_list = [self._unflatten_pair(Y_gbar[a]) for a in range(Y_gbar.shape[0])]

        g_arr = np.stack(g_list, axis=0) if g_list else np.zeros((0, M, M))
        gbar_arr = np.stack(gbar_list, axis=0) if gbar_list else np.zeros((0, M, M))

        # ------------------------------------------------------------------
        # Factor D and Q into lowercase d and q coefficients.
        #
        # Each row of Yd/Yq gives one O_D or O_Q generator.
        # ------------------------------------------------------------------
        if self.spin_explicit:
            Yd_up, rD_up = self._psd_factor(D_up)
            Yd_dn, rD_dn = self._psd_factor(D_dn)
            Yq_up, rQ_up = self._psd_factor(Q_up)
            Yq_dn, rQ_dn = self._psd_factor(Q_dn)

        else:
            Yd, rD = self._psd_factor(D)
            Yq, rQ = self._psd_factor(Q)

        coeffs = dict(
            g=g_arr,
            gbar=gbar_arr,
            d_up=Yd_up if self.spin_explicit else Yd,
            d_dn=Yd_dn if self.spin_explicit else np.zeros((0, M)),
            q_up=Yq_up if self.spin_explicit else Yq,
            q_dn=Yq_dn if self.spin_explicit else np.zeros((0, M)),
        )

        # ------------------------------------------------------------------
        # Sparse support estimate from lowercase generator coefficients.
        # ------------------------------------------------------------------
        sos_support = self._estimate_sparse_d_from_generator_coeffs(
            coeffs,
            threshold=self.support_threshold,
        )

        # ------------------------------------------------------------------
        # Residual check for G11.
        # ------------------------------------------------------------------
        PP = self._diag_pair_indices()

        def term_2(i, j):
            rij = self.pidx(j, i)
            i_j = self.pidx(i, j)

            return 2.0 * (
                np.sum(Gp[rij, PP])
                + np.sum(Gpp[PP, rij])
                - np.sum(Gppp[PP, i_j])
                - np.sum(Gppp[i_j, PP])
            )

        if self.spin_explicit:
            res_up = np.empty((M, M))
            res_dn = np.empty((M, M))

            for i in range(M):
                for j in range(M):
                    lhs_up = term_2(i, j) + (D_up[j, i] - Q_up[j, i])
                    lhs_dn = term_2(i, j) + (D_dn[j, i] - Q_dn[j, i])
                    res_up[i, j] = lhs_up - self.h1[i, j]
                    res_dn[i, j] = lhs_dn - self.h1[i, j]

            res_one_e = max(
                float(np.max(np.abs(res_up))),
                float(np.max(np.abs(res_dn))),
            )

        else:
            res = np.empty((M, M))

            for i in range(M):
                for j in range(M):
                    lhs = term_2(i, j) + (D[i, j] - Q[j, i])
                    res[i, j] = lhs - self.h1[i, j]

            res_one_e = float(np.max(np.abs(res)))

        # ------------------------------------------------------------------
        # Diagnostics and rank accounting.
        # ------------------------------------------------------------------
        diags = dict(
            status=status,
            message=None,
            rank_S=rank_S,
            residual_one_e_maxabs=float(res_one_e),
            g10_rows_dropped=int(self._last_num_g10_dropped),
        )

        if self.spin_explicit:
            diags.update(
                rank_D_up=rD_up,
                rank_D_dn=rD_dn,
                rank_Q_up=rQ_up,
                rank_Q_dn=rQ_dn,
            )
            total_rank = rank_S + rD_up + rD_dn + rQ_up + rQ_dn

        else:
            diags.update(
                rank_D_up=rD,
                rank_D_dn=0,
                rank_Q_up=rQ,
                rank_Q_dn=0,
            )
            total_rank = rank_S + rD + rQ

        result = dict(
            ESOS=float(ESOS),
            lambda_sqrt2=float(lambda_sqrt2),
            Lambda=float(Lambda_eq9),
            coeffs=coeffs,
            diagnostics=diags,
            total_rank=int(total_rank),
            sos_support=sos_support,
        )
        self._last_solution = result
        return result

    def configure_resource_estimates(
        self,
        *,
        E_0,
        eps_qpe,
        Lambda=None,
        E_SOS=None,
        M_spin_orbitals=None,
        R=None,
        bits_precision=20,
        k_r=None,
        optimum=True,
        d_sparse=None,
        reduced_diagonal_doublets=True,
        real_coefficients=True,
        include_a16_overheads=True,
        include_qpe_qubits=True,
        use_last_solution=True,
    ):
        """
        Configure resource-estimate parameters on this same object.

        If ``use_last_solution=True`` and ``solve()`` has already succeeded,
        ``Lambda``, ``E_SOS``, ``R``, and ``d_sparse`` are pulled from the SDP
        result unless explicitly supplied.

        Parameters
        ----------
        E_0 : float
            Reference/target ground-state energy used in the QPE estimate.
        eps_qpe : float
            QPE precision target.
        Lambda : float, optional
            Block-encoding normalization. Defaults to ``solve()["Lambda"]``.
        E_SOS : float, optional
            SOS energy. Defaults to ``solve()["ESOS"]``.
        M_spin_orbitals : int, optional
            Number of spin-orbitals. Defaults to ``2 * norb``.
        R : int, optional
            Total SOS rank/factor count. Defaults to ``solve()["total_rank"]``.
        d_sparse : int, optional
            Sparse support. Defaults to ``solve()["sos_support"]["d_sparse"]``.
        """
        sol = self._last_solution if use_last_solution else None

        if sol is not None and sol.get("diagnostics", {}).get("status") not in (
            "optimal",
            "optimal_inaccurate",
        ):
            sol = None

        if Lambda is None and sol is not None:
            Lambda = sol.get("Lambda")
        if E_SOS is None and sol is not None:
            E_SOS = sol.get("ESOS")
        if R is None and sol is not None:
            R = sol.get("total_rank")
        if d_sparse is None and sol is not None and sol.get("sos_support") is not None:
            d_sparse = sol["sos_support"].get("d_sparse")

        if M_spin_orbitals is None:
            M_spin_orbitals = 2 * self.M

        if Lambda is None:
            raise ValueError("Lambda must be provided or available from a successful solve().")
        if E_SOS is None:
            raise ValueError("E_SOS must be provided or available from a successful solve().")
        if R is None:
            raise ValueError("R must be provided or available from a successful solve().")

        self.resource_Lambda = float(Lambda)
        self.resource_E_SOS = float(E_SOS)
        self.resource_E_0 = float(E_0)
        self.resource_eps_qpe = float(eps_qpe)

        self.M_spin_orbitals = int(M_spin_orbitals)
        self.resource_R = int(R)
        self.bits_precision = int(bits_precision)
        self.optimum = bool(optimum)

        self.d_sparse = None if d_sparse is None else int(d_sparse)

        if k_r is None:
            if self.optimum:
                self.k_r = 0
            else:
                raise ValueError("k_r must be provided when optimum=False.")
        else:
            self.k_r = int(k_r)

        self.reduced_diagonal_doublets = bool(reduced_diagonal_doublets)
        self.real_coefficients = bool(real_coefficients)
        self.include_a16_overheads = bool(include_a16_overheads)
        self.include_qpe_qubits = bool(include_qpe_qubits)

        if self.M_spin_orbitals <= 0:
            raise ValueError("M_spin_orbitals must be positive.")
        if self.resource_R <= 0:
            raise ValueError("R must be positive.")
        if self.M_spin_orbitals % 2 != 0:
            raise ValueError("M_spin_orbitals must be even.")
        if self.k_r < 0:
            raise ValueError("k_r must be a nonnegative integer.")
        if self.bits_precision < 0:
            raise ValueError("bits_precision must be nonnegative.")
        if self.resource_eps_qpe <= 0:
            raise ValueError("eps_qpe must be positive.")
        if self.d_sparse is not None and self.d_sparse <= 0:
            raise ValueError("d_sparse must be positive when provided.")

        self.resource_n_spatial = self.M_spin_orbitals // 2

        if self.d_sparse is not None:
            self.resource_L = int(self.d_sparse)
        else:
            self.resource_L = self._terms_per_O_alpha()

        self._resource_configured = True
        return self

    def estimate_resources(self, **kwargs) -> dict:
        """
        Convenience wrapper: configure resource parameters and return summary().

        Example
        -------
        >>> est = SOSSA_resource_estimates(norb, h1, h2)
        >>> sdp = est.solve()
        >>> resources = est.SDP_resources(E_0=-75.0, eps_qpe=1e-3)
        """
        self.configure_resource_estimates(**kwargs)
        return self.summary()

    def _require_resource_config(self):
        if not getattr(self, "_resource_configured", False):
            raise ValueError(
                "Resource estimates are not configured. Call "
                "configure_resource_estimates(...) or estimate_resources(...)."
            )

    @staticmethod
    def _ceil_log2(x: int) -> int:
        x = int(x)
        if x <= 1:
            return 0
        return int(np.ceil(np.log2(x)))

    @staticmethod
    def _two_adic_valuation(x: int) -> int:
        x = int(x)
        if x <= 0:
            return 0

        eta = 0
        while x % 2 == 0:
            eta += 1
            x //= 2
        return eta

    def _terms_per_O_alpha(self) -> int:
        """
        Dense fallback support count for one Eq.-22 O_alpha.

        Only used when d_sparse is not supplied. Here M is the number of
        spin-orbitals, kept separate from self.M, the SDP spatial-orbital count.
        """
        self._require_resource_config()
        M = self.M_spin_orbitals
        n = self.resource_n_spatial

        if self.reduced_diagonal_doublets:
            return int(1 + 2 * M + M + 8 * n * (n - 1))

        return int(1 + 2 * M + 8 * n * n)

    # ------------------------------------------------------------------
    # Table sizes
    # ------------------------------------------------------------------

    def oalpha_table_size(self) -> int:
        self._require_resource_config()
        return int(self.resource_L)

    def hsqrt_table_size(self) -> int:
        self._require_resource_config()
        if self.d_sparse is not None:
            return int(self.d_sparse)
        return int(self.resource_R * self.resource_L)

    def table_size(self) -> int:
        return self.hsqrt_table_size()

    # ------------------------------------------------------------------
    # Descriptor / QROM output widths
    # ------------------------------------------------------------------

    def local_term_descriptor_bits(self) -> int:
        self._require_resource_config()
        return int(2 * self._ceil_log2(max(2, self.resource_n_spatial)) + 5)

    def hsqrt_term_descriptor_bits(self) -> int:
        self._require_resource_config()
        return int(
            self._ceil_log2(max(2, self.resource_R))
            + self.local_term_descriptor_bits()
        )

    def qrom_output_bits(self, *, include_alpha: bool = True) -> int:
        self._require_resource_config()
        if include_alpha:
            desc = self.hsqrt_term_descriptor_bits()
        else:
            desc = self.local_term_descriptor_bits()

        if self.real_coefficients:
            phase_bits = 1
        else:
            phase_bits = self.bits_precision

        return int(self.bits_precision + 2 * desc + 2 * phase_bits)

    # ------------------------------------------------------------------
    # Effective query count
    # ------------------------------------------------------------------

    def lambda_eff(self) -> float:
        self._require_resource_config()
        E_gap = abs(abs(self.resource_E_0) - self.resource_E_SOS)
        return float(np.sqrt(E_gap * abs(2.0 * self.resource_Lambda - E_gap)))

    def num_walk_steps(self) -> int:
        self._require_resource_config()
        return int(
            np.ceil(np.pi * self.lambda_eff() / (2.0 * self.resource_eps_qpe))
        )

    # ------------------------------------------------------------------
    # QROAM primitives
    # ------------------------------------------------------------------

    @staticmethod
    def _qroam_forward(table_size: int, b_out: int, k: int) -> int:
        fanout = 2 ** int(k)
        return int(np.ceil(table_size / fanout) + (fanout - 1) * b_out)

    @staticmethod
    def _qroam_inverse(table_size: int, k: int) -> int:
        fanout = 2 ** int(k)
        return int(np.ceil(table_size / fanout) + fanout)

    def _qroam_pair_toffoli(self, table_size: int, b_out: int, k: int) -> int:
        return int(
            self._qroam_forward(table_size, b_out, k)
            + self._qroam_inverse(table_size, k)
        )

    def _choose_k_opt(self, table_size: int, b_out: int) -> int:
        table_size = int(table_size)
        b_out = int(b_out)

        if table_size <= 1:
            return 0

        k_max = int(np.floor(np.log2(table_size)))

        k0 = 0.5 * np.log2(max(1.0, 2.0 * table_size / max(1.0, b_out + 1)))
        k0 = int(np.clip(np.round(k0), 0, k_max))

        cand = sorted(
            set([0, k_max, k0 - 3, k0 - 2, k0 - 1, k0, k0 + 1, k0 + 2, k0 + 3])
        )
        cand = [k for k in cand if 0 <= k <= k_max]

        best_k = cand[0]
        best_val = self._qroam_pair_toffoli(table_size, b_out, best_k)

        for k in cand[1:]:
            val = self._qroam_pair_toffoli(table_size, b_out, k)
            if val < best_val:
                best_val = val
                best_k = k

        return int(best_k)

    def chosen_k_for_table(self, table_size: int, b_out: int) -> int:
        self._require_resource_config()
        if not self.optimum:
            return self.k_r
        return self._choose_k_opt(table_size, b_out)

    def chosen_k(self) -> int:
        return self.chosen_k_for_table(
            self.hsqrt_table_size(),
            self.qrom_output_bits(include_alpha=True),
        )

    # ------------------------------------------------------------------
    # SELECT and block-encoding costs
    # ------------------------------------------------------------------

    def select_toffoli(self) -> int:
        self._require_resource_config()
        return int(max(0, 2 * self.M_spin_orbitals - 2))

    def _block_encoding_toffoli_for_table(
        self,
        *,
        table_size: int,
        b_out: int,
    ) -> int:
        self._require_resource_config()
        table_size = int(table_size)
        b_out = int(b_out)

        k = self.chosen_k_for_table(table_size, b_out)

        qroam = self._qroam_pair_toffoli(table_size, b_out, k)
        sel = self.select_toffoli()

        if not self.include_a16_overheads:
            return int(qroam + sel)

        n_d = self._ceil_log2(max(2, table_size))
        eta = self._two_adic_valuation(table_size)
        b_r = self.bits_precision

        overhead = (
            2 * self.bits_precision
            + 7 * n_d
            - 6 * eta
            + 4 * b_r
            - 19
        )

        return int(qroam + sel + max(0, overhead))

    def hsqrt_block_encoding_toffoli(self) -> int:
        return self._block_encoding_toffoli_for_table(
            table_size=self.hsqrt_table_size(),
            b_out=self.qrom_output_bits(include_alpha=True),
        )

    def block_encoding_toffoli(self) -> int:
        return self.hsqrt_block_encoding_toffoli()

    # ------------------------------------------------------------------
    # Walk and total Toffoli costs
    # ------------------------------------------------------------------

    def walk_reflection_toffoli(self) -> int:
        self._require_resource_config()
        T = self.hsqrt_table_size()
        n_d = self._ceil_log2(max(2, T))
        return int(n_d + self.bits_precision + 2)

    def walk_toffoli(self) -> int:
        return int(
            2 * self.hsqrt_block_encoding_toffoli()
            + self.walk_reflection_toffoli()
        )

    def toffoli_count(self) -> int:
        return int(self.num_walk_steps() * self.walk_toffoli())

    # ------------------------------------------------------------------
    # Qubit counts
    # ------------------------------------------------------------------

    def _qubit_count_for_table(
        self,
        *,
        table_size: int,
        b_out: int,
        include_qpe: bool,
    ) -> int:
        self._require_resource_config()
        table_size = int(table_size)
        b_out = int(b_out)

        k = self.chosen_k_for_table(table_size, b_out)
        fanout = 2 ** k

        qrom_workspace = int(
            b_out * fanout
            + self._ceil_log2(max(2, int(np.ceil(table_size / fanout))))
        )

        prepare_address = self._ceil_log2(max(2, table_size))

        if include_qpe:
            I = self.num_walk_steps()
            qpe_qubits = 2 * self._ceil_log2(max(2, I + 1)) - 1
        else:
            qpe_qubits = 0

        return int(
            self.M_spin_orbitals
            + qpe_qubits
            + prepare_address
            + self.bits_precision
            + self.bits_precision
            + qrom_workspace
            + 3
        )

    def oalpha_qubit_count(self) -> int:
        return self._qubit_count_for_table(
            table_size=self.oalpha_table_size(),
            b_out=self.qrom_output_bits(include_alpha=False),
            include_qpe=False,
        )

    def qubit_count(self) -> int:
        return self._qubit_count_for_table(
            table_size=self.hsqrt_table_size(),
            b_out=self.qrom_output_bits(include_alpha=True),
            include_qpe=self.include_qpe_qubits,
        )

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------

    def support_summary(self) -> dict:
        self._require_resource_config()
        return {
            "using_generator_coefficient_d_sparse": self.d_sparse is not None,
            "d_sparse": self.d_sparse,
            "dense_fallback_L": None if self.d_sparse is not None else self.resource_L,
            "R": self.resource_R,
            "M_spin_orbitals": self.M_spin_orbitals,
            "n_spatial_orbitals": self.resource_n_spatial,
            "table_size": self.hsqrt_table_size(),
            "chosen_k": self.chosen_k(),
            "qrom_output_bits": self.qrom_output_bits(include_alpha=True),
        }

    def summary(self) -> dict:
        self._require_resource_config()
        return {
            **self.support_summary(),
            "Lambda": self.resource_Lambda,
            "E_SOS": self.resource_E_SOS,
            "E_0": self.resource_E_0,
            "eps_qpe": self.resource_eps_qpe,
            "lambda_eff": self.lambda_eff(),
            "num_walk_steps": self.num_walk_steps(),
            "block_encoding_toffoli": self.block_encoding_toffoli(),
            "walk_reflection_toffoli": self.walk_reflection_toffoli(),
            "walk_toffoli": self.walk_toffoli(),
            "total_toffoli": self.toffoli_count(),
            "logical_qubits": self.qubit_count(),
        }

