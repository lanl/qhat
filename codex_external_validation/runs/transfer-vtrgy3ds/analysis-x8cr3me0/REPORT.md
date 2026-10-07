# Prospective transfer results

Campaign: `transfer-vtrgy3ds`.
Four new molecular families; two idealized geometries per family; 6-31g/CAS(4e,5o), HF.
Rows/geometries/times are not independent molecules. No chemically converged molecular accuracy claim.

Completed 216 baseline + 56 control rows. Prediction floor exclusions: 1.
Primary 5% prediction failures: 0/143.
Max relative I error r=100: 0.101498%; r=200: 0.050571%.
Rank agreement (>1% actual margin, r=100): dynamic 72/72, static HF-BCH 72/72; excluded pairs 0.
Median relative error (r=100,200): dynamic 0.002302%, static 4.337%.

## All T=1,r=100 anchors (E=1-|overlap|)

| Case | Fermionic signed | JW signed | JW magnitude | Best |
|---|---:|---:|---:|---|
| HCN_s1.00 | 3.353726552e-07 | 8.033052188e-08 | 2.949983663e-07 | jw_signed |
| HCN_s1.25 | 1.006384443e-07 | 4.081940468e-07 | 1.105762449e-07 | fermionic_signed |
| CH2O_s1.00 | 1.871830644e-07 | 1.576213833e-07 | 4.814401578e-07 | jw_signed |
| CH2O_s1.25 | 1.277840939e-07 | 2.417213055e-07 | 3.619261166e-07 | fermionic_signed |
| H2S_s1.00 | 4.318094942e-08 | 4.678993427e-09 | 4.262984733e-08 | jw_signed |
| H2S_s1.25 | 1.241583806e-07 | 1.941413440e-08 | 8.412056041e-08 | jw_signed |
| SiH4_s1.00 | 4.406921782e-08 | 3.452008887e-10 | 4.926888003e-08 | jw_signed |
| SiH4_s1.25 | 6.668182130e-08 | 5.219538379e-09 | 6.808574889e-08 | jw_signed |

## Interpretation guards

Predictions are unfitted leading-error propagation, a standard perturbative construction, not a new general theory.
The static comparator is a local HF BCH proxy, not an implementation of published best algorithms.
Symmetry preservation is not a guarantee of lower total state error. Retain every contrary ranking.
Baseline C range at r=100: [83.19762988625874, 99.77174352675672]; C<1 cases: 0. A C reduction above 1 is NOT net cancellation.

## Machine-readable checks

```json
{
  "campaign": "/Users/albertlee0125/Repos/qhat/codex_external_validation/runs/transfer-vtrgy3ds",
  "unique_families": 4,
  "geometries": 8,
  "baseline_rows": 216,
  "control_rows": 56,
  "prediction_rows": 216,
  "excluded_below_floor_rows": 1,
  "primary_eligible_rows_r100_r200": 143,
  "primary_prediction_failures": 0,
  "r100_dynamic_max_relative_error": 0.00101497557320962,
  "r200_dynamic_max_relative_error": 0.0005057146776232813,
  "dynamic_median_relative_error": 2.302050295244218e-05,
  "static_median_relative_error": 0.04337483607972126,
  "static_max_relative_error": 0.2144040656210584,
  "rank_eligible_pairs": 72,
  "rank_excluded_pairs": 0,
  "dynamic_rank_correct": 72,
  "static_rank_correct": 72,
  "case_time_exact_winner_counts": {
    "fermionic_signed": 6,
    "jw_signed": 18,
    "jw_magnitude": 0
  },
  "r_scaling_slope_min_max": [
    -2.0006298067899486,
    -1.9988921858678717
  ],
  "max_norm_drift": 1.311839525897085e-12,
  "max_defect_reconstruction": 5.144140275100358e-14,
  "max_overlap_reconstruction": 5.0104627708058945e-14,
  "max_projected_reconstruction_difference": 1.107717924683652e-17,
  "fermionic_max_spin_leakage": 5.04172425712769e-37,
  "jw_signed_max_spin_leakage": 1.172947950545761e-07,
  "baseline_C_min_max": [
    83.19762988625874,
    99.77174352675672
  ],
  "baseline_net_destructive_count": 0,
  "source_sha256": "e0e072fcb598cd5b15b7fdf8161d1b0d69c3a7290190c7c0d483b986f4adcf23",
  "repeat": {
    "path": "/Users/albertlee0125/Repos/qhat/codex_external_validation/runs/transfer-hux1k985",
    "all_272_result_rows_exactly_equal": false,
    "all_physical_metrics_within_rtol1e_8_atol1e_20": true,
    "physical_metric_comparisons": 2288,
    "max_relative_difference_above_1e_20": 3.2221079599655757e-13,
    "all_216_predictions_exactly_equal": true,
    "all_input_and_order_hashes_equal": true
  }
}
```
