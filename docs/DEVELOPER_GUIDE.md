# QHAT developer's guide

For people changing QHAT's code, especially its Qualtran bloqs. User-facing options are in
[`analysis/README.md`](../analysis/README.md); this guide covers how the code is built and what
goes wrong when it is changed carelessly.

Contents: [Writing a bloq](#1-writing-a-bloq) · [Pitfalls](#2-pitfalls) ·
[Cost-model conventions](#3-cost-model-conventions) · [Testing](#4-testing-conventions) ·
[Environment](#5-environment) · [AI usage policy](#6-ai-usage-policy)

## 1. Writing a bloq

QHAT's own bloqs are `Trotterization` (`common/trotter_flattened.py`),
`CommutingPauliStringEvolution`, `PauliStringEvolution`, and `FlexibleQPE`. Controls live in a
shared helper module, `common/controls.py`; read its docstring first.

The idea behind all of them: Qualtran's generic fallbacks (wrapping in `Controlled` or `Adjoint`,
decomposing to find a count) are correct but slow, and they lose our optimizations. So each bloq
implements the pieces below itself. Qualtran's default is fine in some cases and silently bad in
others; the table says which.

| Piece | What it does | Qualtran's default | Our rule |
|---|---|---|---|
| `signature` | Registers. Controls are a register `ctrl` of shape `(len(ctrl_values),)` placed before `q`. | none (required) | Implement. |
| `build_composite_bloq` | The exact decomposition. | none (required, or implement `decompose_bloq`) | Must equal the matrix **exactly**, global phase included, because a global phase becomes a relative phase under control (see pitfalls). Build rotations from Qualtran's `Rz`, not Cirq's `Z**x`. |
| `_t_complexity_` | Fast T, Clifford and rotation counts for `t_complexity` and pyLIQTR. | decomposes the bloq | Override for any bloq with a long or repetitive decomposition. It must agree with the decomposition (tests check this). |
| `build_call_graph` | The callees and their counts, used by `QECGatesCost` and call-graph tools. | decomposes the bloq | Override for long decompositions. Return the **same** bloqs with the same parameters as the decomposition, never a callee with count 0. Check with `dict(b.build_call_graph(None)) == dict(b.decompose_bloq().build_call_graph(None))`. |
| `my_static_costs` | A cost without decomposing. We use it for `QubitCount`. | decomposes | Override when the peak is known (system qubits + controls + ancillas). |
| `add_my_tensors` / `tensor_contract` | The unitary, from an exact matrix. | contracts the decomposition gate by gate | Override: it is much faster and is the reference the decomposition is tested against. |
| `adjoint` | U†. | `Adjoint(bloq)` wrapper | Flip an explicit `is_adjoint` field. Decomposition and tensors honour it. |
| `get_ctrl_system` | What `.controlled()` returns. | `Controlled(bloq)` wrapper | Return our own bloq with a new value prepended to `ctrl_values` (see `controls.get_ctrl_system`). |
| `__pow__` | `U**k`. | none | Positive k as before; `U**(-k)` is `(U**k).adjoint()`; 0 and non-integers raise. |

Details worth knowing:

- **Controls are fields.** `ctrl_values=(1, 0)` means "act when the first control is |1> and the
  second is |0>". `.controlled()` prepends one value, so chains such as
  `U.controlled().adjoint().controlled()` stay one bloq of our type. Only single-qubit controls
  on |1> or |0> are supported; other `CtrlSpec`s raise.
- **Adjoint.** Don't negate `time` to get U†. That is correct only for symmetric (palindromic)
  Trotter sequences; first and third order are not. `U(-t)` and `U(t).adjoint()` are deliberately
  different objects. Because `time` can therefore be negative, any code that reads it as a
  duration (precision, step size, error bounds) must use `abs(time)`.
- **Equality and hashing.** Every field that changes the circuit (`is_adjoint`, `ctrl_values`,
  `multi_control_method`, `controlled_rotation_method`) is part of `eq`/hash.
- **Where the control goes.** In a Pauli-string evolution only the central rotation needs the
  control: the basis changes and the CNOT ladder cancel when it is skipped. A control on |0> needs
  no X gates; it flips the sign of the second half-rotation.
- **Several controls.** AND them into an ancilla once, at the outermost bloq (`"ancilla"`
  method); inner bloqs then see one control. `"direct"` uses no AND ancilla and makes each
  rotation multi-controlled instead (2 rotations plus relative-phase Toffolis, borrowing idle
  qubits as dirty ancillas, or a 2^k Gray-code walk when too few are idle).
- **Controlled rotation.** The default (`"textbook"`) is Rz(θ/2) · CNOT · Rz(−θ/2) · CNOT.
  `"ancilla"` rotates an AND ancilla instead (1 rotation + 4 T). It leaves a phase on the control,
  which commutes with everything controlled on that qubit, so a parent bloq sums its terms'
  phases and applies them once (`omit_control_phase`, `controls.add_controlled_phase`). Without
  that merge the option costs more than the textbook circuit.
- **Rotation counting.** Every Pauli term counts as one rotation whatever its angle; zero
  coefficients and multiples of π/4 get no special treatment. An all-identity string is a
  `GlobalPhase` (0 rotations uncontrolled, 1 under control).
- **Empty decompositions.** In Qualtran 0.4, `t_complexity` of an empty `CompositeBloq`
  raises. A bloq whose decomposition can be empty must count itself.

### Adding an option to the Trotter bloqs

Put the field on all three bloqs and in `eq`/hash; validate it in `controls.py`; add it to
`UnitaryConfiguration` (`analysis/config_types.py`) and pass it on in `analysis/unitary.py`;
reject it with `trotter_implementation = "original"`; and record it in the TOML only if the user
set it (see section 3).

## 2. Pitfalls

Each of these has bitten this code before.

**Phases**
- **Cirq's `Z**x` is not Rz.** It equals Rz(θ) times a global phase, and under control that phase
  becomes a relative phase. A controlled decomposition built from it applies e^{iφ}·U, not U.
  Matrices from `tensor_contract` stay correct, so only decompositions (and anything that uses
  them, like pyLIQTR) are wrong. Use Qualtran's `Rz` (exactly exp(−iθZ/2)).
- **Never build a Pauli-string gate from `cirq.PauliStringPhasorGate`.** Cirq 1.4 treats `I` as
  `Z`. Act only on the non-identity qubits.
- **`SGate(is_adjoint=True)` converts to plain `cirq.S` in Qualtran 0.4**, so any string with a Y
  came out wrong in Cirq and pyLIQTR. Use `CirqSafeSGate` (`common/pauli_string_evolution.py`) in
  any decomposition that may go through Cirq.

**Wrappers that lose our optimizations**
- **Add controls with `.controlled()`, never `Controlled(...)`.** Qualtran's `Controlled` around
  one of our bloqs controls each sub-bloq separately. That repeats the AND ladder for every
  term, and it loses the merged control phase of the `"ancilla"` rotation method (one correction
  per term instead of one per bloq). The unitary stays exact, so only the counts reveal it.
- **`TextbookQPE` always uses `Controlled`.** It builds `Controlled(unitary)` directly, so it never
  reaches our `get_ctrl_system`. Use `FlexibleQPE`. For the same reason
  `controlled_rotation_method = "ancilla"` with `"QPE: qualtran textbook"` raises an error.
  Any new path that wraps our bloqs in `Controlled` must reject settings that depend on the merge.
- **The default call graph is slow.** It builds the whole flat circuit, so it grows with the
  number of steps (about 25 s for `QECGatesCost` of a 500-step QPE). Override `build_call_graph`
  and `_t_complexity_`, and speed up the leaf bloqs too; with many distinct terms their
  decompositions dominate.
- **pyLIQTR uses `_t_complexity_` when a bloq has one**, and decomposes the bloqs that don't
  (like `Controlled`). Cost changes to decompositions can therefore move pyLIQTR numbers even when
  `_t_complexity_` is unchanged. pyLIQTR's `LogicalQubits` counts registers only and is absent
  for some bloqs; QHAT's `qubit_count` always comes from `QubitCount`.

**Measuring things**
- **`get_cost_value(composite_bloq, QubitCount())` returns the register width**, not the peak
  of its contents. To check a bloq's `my_static_costs` against its decomposition, use
  `_cbloq_max_width` (see `decomposition_width` in `common/tests/test_multi_control.py`).
- **`QECGatesCost` and `t_complexity` disagree** on `And` (one `and_bloq` versus 4 T + 9
  Cliffords; the uncompute is 1 Clifford + a measurement versus 4 Cliffords) and on `CZPowGate`,
  which `t_complexity` counts as a rotation and `QECGatesCost` doesn't. Qualtran's textbook
  inverse QFT uses `CZPowGate`, so the two estimators differ by `n(n−1)/2` rotations for `n` phase
  qubits. QHAT takes rotations from `t_complexity`.
- **`CompositeBloq.bloq_instances` is an unordered set.** To check gate order, use
  `iter_bloqnections()`.
- **`cirq.unitary(bloq)` raises.** A bloq isn't a Cirq gate. Use `BloqAsCirqGate(bloq)` or
  `bloq.tensor_contract()`. `cirq.unitary` also can't contract the `"ancilla"` multi-control
  circuit (it has ancillas), nor the measurement-based `And†`.
- **pyLIQTR's rotation conversion calls `random.seed(0)`**, resetting Python's global `random`
  state. Nothing in QHAT relies on that state; don't start.

**Known issue.** The Qualtran estimator crashes on `"QPE: qualtran textbook"`:
`get_cost_value(TextbookQPE(...), QubitCount())` raises an `AssertionError` in Qualtran's Cirq
interop, because `TextbookQPE` has no decomposition and falls back to Cirq. pyLIQTR's estimator
works, and so does `"QPE: QHAT flexible"`. Nothing fixes it yet.

## 3. Cost-model conventions

- **Qubit count means the peak.** `qubit_count` is the largest number of qubits in use at any
  point, including ancillas allocated and freed inside the circuit, not the register width.
  `QubitCount` is the source of truth. The log prints a qubit summary (`analysis/qubit_summary.py`);
  if its components don't add up to the total, it warns.
- **Rotations are totalled, then converted once.** Count each rotation as one rotation, add them
  over the whole circuit, and convert the total to T + Clifford at the end with
  `convert_rotations` (`analysis/resource_estimation.py`). Never convert per bloq. The model is
  pyLIQTR's (3.02·log2(1/ε) + 0.77 T per rotation, 2 Cliffords per T + 1, with a seeded random
  offset), the same for both estimators, so the T count is reproducible but not smooth in ε.
  The precision is `analysis.rotation_gate_precision` (per rotation, default 1e-10) or
  `analysis.algorithm_precision` (divided by the total rotation count), never both. Results carry
  both the converted counts and `T_count_unconverted` / `Clifford_count_unconverted`.
- **Number of controls.** `algorithm.control_values` (e.g. `[1, 0, 1]`) sets the controls of
  `"controlled time evolution"`; other methods raise if it is set. Unset means one control on |1>.
- **TOML versus log.** The TOML records what the user set; the log records what QHAT chose. An
  option the user did not set must not appear in the TOML. The log says which multi-control
  method was used, how many ancillas it added, and when an option is ignored (for example
  `multi_control_method` with fewer than 2 controls).
- **Conflicting settings raise errors**, and so do settings that leave nothing to choose from.
  Honour every option the user set. Validate right after reading the configuration, before any
  expensive work; `check_algorithm_settings` (`analysis/algorithm.py`) is where algorithm-level
  conflicts go, and `build_algorithm` calls it so every caller is covered. Examples:
  `trotter_implementation = "original"` with any controlled algorithm; `controlled_power_method =
  "sign"` with an asymmetric Trotter or with `controlled_rotation_method = "ancilla"`;
  `controlled_rotation_method = "ancilla"` with `multi_control_method = "direct"` and 2 or more
  controls.
- **Trotter-only options with other unitaries** (Pauli LCU, double factorization) may be silently
  ignored. When a Trotterization is used and the option affects it, the TOML must record it.
- **Which options exist.** `multi_control_method` (`"ancilla"` default / `"direct"`),
  `controlled_rotation_method` (`"textbook"` default / `"ancilla"`) and `algorithm.
  controlled_power_method` (`"controlled"` default / `"sign"`, "QPE: QHAT flexible" only). The
  defaults are chosen to be the safe, cheapest-in-qubits option; users opt into the others.

## 4. Testing conventions

- **Check exactness against explicitly built matrices.** A controlled unitary is U on the block
  matching the control pattern and the identity elsewhere (controls most significant), compared
  exactly, not "equal up to a phase". Phase-blind tests miss exactly the bugs in the pitfalls
  above. See `block_matrix` and `check_exact` in `common/tests/test_multi_control.py`, which
  checks both `tensor_contract()` and the decomposition.
- **Use non-special parameters.** Times that aren't 2π, more than one Trotter step, several
  Trotter orders (including non-symmetric ones such as first and third), and terms that don't
  commute. Commuting terms hide reversal bugs; a symmetric sequence hides adjoint bugs.
- **Prefer real input files** (`analysis/tests/*.hdf5`, `*.fci`, `xxz_cs2cocl4_10.json`,
  `julia_trotter/Li-Li_*`) where practical.
- **Hand-count the costs.** Tests state the expected T, Clifford and rotation counts for small
  strings, and check that `_t_complexity_`, `build_call_graph` and the decomposition agree.
- **Compare numbers before and after** any change to a decomposition. Counts that pass through a
  generic wrapper depend strongly on how the wrapped bloq decomposes.
- **Run the tests** with `python3.11 -m pytest common/tests analysis/tests` (about 60 s; over
  2,700 tests). Run the files you touched while working and the whole suite before finishing.
  Anything that may take over 10 minutes (impact tables over many inputs) goes in the background.
- **CI** runs on `ubuntu-24.04`, pinned in `.github/workflows/python-app.yml`, because pyLIQTR
  needs Julia 1.9, which doesn't run on newer runner images. If only the `analysis-pauli-lcu` and
  `analysis-dbl-factor` jobs fail, with a `juliapkg` traceback, check the runner image before
  suspecting the code. GitHub's support for `ubuntu-24.04` is expected to end around April 2027.

## 5. Environment

- Use `python3.11` for everything, and a virtual environment per worktree
  (`python3.11 -m venv venv && venv/bin/pip install -e .`). An editable install points at one
  directory, so a shared environment makes a worktree silently import another checkout's code.
  Check with `python -c "import qhat; print(qhat.__path__)"`.
- Stay on Qualtran 0.4.0, Cirq 1.4.0 and pyLIQTR 1.4.2: pyLIQTR blocks a newer Qualtran. An
  upgrade is a separate piece of work. Keep our overrides small and in one place per bloq so it
  is cheap to revisit.

## 6. AI usage policy

AI tools are allowed. The human contributor is responsible for verifying AI-assisted work (code,
tests, analysis, citations, resource numbers) as if they had written it.

- Citations are copied from the source itself, never from an AI's memory.
- Merge requests disclose AI use.
- Name things descriptively (`num_phase_qubits`), not with one paper's symbols.
- The team sets this policy; no approval beyond the team is needed.
