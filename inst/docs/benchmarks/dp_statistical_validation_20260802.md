# dsVert DP statistical validation — 2026-08-02

## Scope

- fixed finite-dataset distributional simulation using the ideal two-peer geometric convolution; not DSI E2E and not a productive cryptographic-sampler replay.
- DP mechanism noise and public grid/quantization only; no sampling or population confidence intervals.
- Fixed mechanism replicates per family: `1000`.
- Seeds: `describe=20260802`, `epidemiology=20260803`, `mantel_haenszel=20260807`, `survival=20260804`, `gaussian=20260805`, `spectral=20260806`.
- Runtime is descriptive only and has no pass/fail threshold.
- The central oracle is the bounded finite dataset: pre-quantization
  moments for fidelity, and fixed public-bin endpoints for quantiles.
  Zero-noise identities separately target the deployed quantized
  estimand, and its reported grid region must contain the bounded truth.
  No population or sampling confidence interval is invented.

## Gates

| gate | pass | detail |
| --- | --- | --- |
| deployed_accuracy_contract_at_least_95_percent | TRUE | minimum=0.950026445836 |
| truth_covered_whenever_certified_coordinate_box_holds | TRUE | checked=63493 |
| all_main_regions_ordered_and_finite | TRUE | rows=70000; width_range=[0,4] |
| all_point_estimates_finite | TRUE | rows=70000; all finite |
| central_oracle_zero_noise_identity | TRUE | max_abs_error=4.44089e-16 |
| gaussian_quantized_oracle_fidelity | TRUE | max_abs_error=0.000121832 |
| correlation_pca_quantized_oracle_fidelity | TRUE | max_abs_error=0.000817716 |
| describe_quantization_region_contains_bounded_truth | TRUE | mean and variance bounded truth inside no-noise quantization region |
| quantile_median_identity | TRUE | max_abs_error=0 |
| epidemiology_diagnostic_and_causal_standardization_identities | TRUE | max_abs_error=1.332268e-15 |
| mantel_haenszel_zero_cost_no_classical_dp_inference | TRUE | oracle_error=0; zero_cost_no_classical=TRUE |
| survival_competing_risk_and_monotonicity_identities | TRUE | max_abs_error=3.330669e-16 |
| gaussian_solve_algebra_identity | TRUE | max_abs_error=0 |
| pca_trace_and_psd_identities | TRUE | max_abs_error=5.77316e-15 |
| all_edge_cases_typed | TRUE | passed=16/16 |
| release_only_postprocessors_have_no_dsi_calls | TRUE | static body audit of release-only public postprocessors |

## Deployed accuracy contracts

| family | epsilon | natural_l1_sensitivity | coordinate_count | scale | simultaneous_radius_95 | implementation_tv_upper_bound | certified_failure_probability_upper | certified_coverage_lower | accuracy_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| describe_quantile_median | 6 | 5 | 12 | 256 | 5.82031 | 0.00000000000000000000000000000157772 | 0.0498262 | 0.950174 | exact ideal two-sided-geometric convolution tail with union bound; two-peer finite-sampler TV deducted; fixed-clamp range applied |
| epidemiology_2x2 | 3 | 2 | 4 | 256 | 3.82422 | 0.00000000000000000000000000000157772 | 0.0497965 | 0.950203 | exact ideal two-sided-geometric convolution tail with union bound; two-peer finite-sampler TV deducted; fixed-clamp range applied |
| mantel_haenszel | 3 | 2 | 12 | 256 | 4.65625 | 0.00000000000000000000000000000157772 | 0.0498001 | 0.950200 | exact ideal two-sided-geometric convolution tail with union bound; two-peer finite-sampler TV deducted; fixed-clamp range applied |
| survival_rmst | 3 | 2 | 19 | 256 | 5.00000 | 0.00000000000000000000000000000157772 | 0.0497849 | 0.950215 | exact ideal two-sided-geometric convolution tail with union bound; two-peer finite-sampler TV deducted; fixed-clamp range applied |
| gaussian_correlation_pca | 6 | 12 | 12 | 256 | 13.96484 | 0.00000000000000000000000000000157772 | 0.0499736 | 0.950026 | exact ideal two-sided-geometric convolution tail with union bound; two-peer finite-sampler TV deducted; fixed-clamp range applied |

The analytic gate combines the exact ideal convolution tail, the union
bound over the released coordinates, and the productive sampler's
published two-peer total-variation allowance. Region coverage is also
required deterministically whenever every exact coordinate lies inside
that certified simultaneous box.

## Bias, RMSE, mechanism coverage, width, and runtime

| family | method | estimand | observations | bias | rmse | mechanism_region_available | mechanism_coverage | conditional_coordinate_coverage | mean_region_width | coordinate_event_rate | estimate_finite_rate | finite_rate | elapsed_seconds_family | replicates |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| correlation | ds.vertCor | x1_x2 | 1000 | -0.00327510809 | 0.06282345 | TRUE | 1.000 | 1 | 1.5996355 | 0.956 | 1 | 1 | 2.910 | 1000 |
| correlation | ds.vertCor | x1_y | 1000 | -0.03534219477 | 0.10044264 | TRUE | 1.000 | 1 | 1.9958801 | 0.956 | 1 | 1 | 2.910 | 1000 |
| correlation | ds.vertCor | x2_y | 1000 | -0.04491823984 | 0.09882192 | TRUE | 1.000 | 1 | 1.9881560 | 0.956 | 1 | 1 | 2.910 | 1000 |
| describe | ds.vertDPDescribe | bounded_mean | 1000 | -0.00114369339 | 0.01943015 | TRUE | 1.000 | 1 | 0.1828675 | 0.948 | 1 | 1 | 11.357 | 1000 |
| describe | ds.vertDPDescribe | bounded_population_variance | 1000 | 0.00448184542 | 0.09867139 | TRUE | 1.000 | 1 | 1.4069843 | 0.948 | 1 | 1 | 11.357 | 1000 |
| describe | ds.vertDPMedian | median_public_bin_endpoint | 1000 | 0.02600000000 | 0.11401754 | TRUE | 1.000 | 1 | 1.0000000 | 0.948 | 1 | 1 | 11.357 | 1000 |
| describe | ds.vertDPQuantile | q25_public_bin_endpoint | 1000 | -0.00150000000 | 0.02738613 | TRUE | 1.000 | 1 | 1.0000000 | 0.948 | 1 | 1 | 11.357 | 1000 |
| describe | ds.vertDPQuantile | q50_public_bin_endpoint | 1000 | 0.02600000000 | 0.11401754 | TRUE | 1.000 | 1 | 1.0000000 | 0.948 | 1 | 1 | 11.357 | 1000 |
| describe | ds.vertDPQuantile | q75_public_bin_endpoint | 1000 | 0.01600000000 | 0.08944272 | TRUE | 1.000 | 1 | 1.0000000 | 0.948 | 1 | 1 | 11.357 | 1000 |
| epidemiology | ds.vertDPCausalStandardization | odds_ratio | 1000 | 0.00105894403 | 0.05234335 | TRUE | 1.000 | 1 | 0.5331044 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPCausalStandardization | risk_control | 1000 | -0.00007930139 | 0.00276644 | TRUE | 0.994 | 1 | 0.0191189 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPCausalStandardization | risk_difference | 1000 | -0.00000991032 | 0.00356937 | TRUE | 1.000 | 1 | 0.0382356 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPCausalStandardization | risk_ratio | 1000 | 0.00072712406 | 0.02978529 | TRUE | 0.999 | 1 | 0.2877187 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPCausalStandardization | risk_treated | 1000 | -0.00008921170 | 0.00229135 | TRUE | 0.999 | 1 | 0.0191167 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDiagnostic2x2 | accuracy | 1000 | -0.00001686809 | 0.00166678 | TRUE | 1.000 | 1 | 0.0191176 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDiagnostic2x2 | balanced_accuracy | 1000 | -0.00000495516 | 0.00178468 | TRUE | 1.000 | 1 | 0.0191178 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDiagnostic2x2 | diagnostic_odds_ratio | 1000 | 0.00105894403 | 0.05234335 | TRUE | 1.000 | 1 | 0.5331044 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDiagnostic2x2 | f1_score | 1000 | -0.00005562002 | 0.00244267 | TRUE | 1.000 | 1 | 0.0238988 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDiagnostic2x2 | lr_negative | 1000 | 0.00004601661 | 0.00384038 | TRUE | 1.000 | 1 | 0.0418231 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDiagnostic2x2 | lr_positive | 1000 | 0.00072712406 | 0.02978529 | TRUE | 0.999 | 1 | 0.2877187 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDiagnostic2x2 | npv | 1000 | -0.00003970445 | 0.00161762 | TRUE | 1.000 | 1 | 0.0136540 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDiagnostic2x2 | ppv | 1000 | 0.00006854809 | 0.00403612 | TRUE | 0.999 | 1 | 0.0318722 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDiagnostic2x2 | prevalence | 1000 | 0.00002768967 | 0.00156992 | TRUE | 1.000 | 1 | 0.0191176 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDiagnostic2x2 | sensitivity | 1000 | -0.00008921170 | 0.00229135 | TRUE | 0.999 | 1 | 0.0191167 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDiagnostic2x2 | specificity | 1000 | 0.00007930139 | 0.00276644 | TRUE | 0.994 | 1 | 0.0191189 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPDirectStandardization | direct_standardized_risk | 1000 | -0.00008524758 | 0.00177588 | TRUE | 1.000 | 1 | 0.0191176 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPEpi2x2 | attributable_fraction_exposed | 1000 | 0.00007082668 | 0.00745378 | TRUE | 0.999 | 1 | 0.0717499 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPEpi2x2 | odds_ratio | 1000 | 0.00105894403 | 0.05234335 | TRUE | 1.000 | 1 | 0.5331044 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPEpi2x2 | population_attributable_fraction | 1000 | 0.00010997604 | 0.00690672 | TRUE | 1.000 | 1 | 0.1063444 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPEpi2x2 | population_risk | 1000 | -0.00007671336 | 0.00178850 | TRUE | 1.000 | 1 | 0.0191176 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPEpi2x2 | risk_difference | 1000 | -0.00000991032 | 0.00356937 | TRUE | 1.000 | 1 | 0.0382356 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPEpi2x2 | risk_exposed | 1000 | -0.00008921170 | 0.00229135 | TRUE | 0.999 | 1 | 0.0191167 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPEpi2x2 | risk_ratio | 1000 | 0.00072712406 | 0.02978529 | TRUE | 0.999 | 1 | 0.2877187 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPEpi2x2 | risk_unexposed | 1000 | -0.00007930139 | 0.00276644 | TRUE | 0.994 | 1 | 0.0191189 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPIndirectStandardization | observed_expected_ratio | 1000 | -0.00036432097 | 0.00803412 | TRUE | 1.000 | 1 | 0.0849657 | 0.959 | 1 | 1 | 8.971 | 1000 |
| epidemiology | ds.vertDPROC | auc | 1000 | -0.00000495516 | 0.00178468 | TRUE | 1.000 | 1 | 0.0880347 | 0.959 | 1 | 1 | 8.971 | 1000 |
| gaussian | ds.vertDPGaussian | (Intercept) | 1000 | 0.00807754540 | 0.04820839 | FALSE | NA | NA | NA | 0.958 | 1 | 1 | 0.399 | 1000 |
| gaussian | ds.vertDPGaussian | x1 | 1000 | -0.00825890318 | 0.05965249 | FALSE | NA | NA | NA | 0.958 | 1 | 1 | 0.399 | 1000 |
| gaussian | ds.vertDPGaussian | x2 | 1000 | -0.00770860320 | 0.05677876 | FALSE | NA | NA | NA | 0.958 | 1 | 1 | 0.399 | 1000 |
| mantel_haenszel | ds.vertDPMantelHaenszel | common_odds_ratio | 1000 | -0.00235298728 | 0.07466458 | TRUE | 1.000 | 1 | 1.6910748 | 0.956 | 1 | 1 | 0.983 | 1000 |
| pca | ds.vertPCA | eigenvalue_PC1 | 1000 | -0.05250690179 | 0.10767117 | TRUE | 1.000 | 1 | 3.9931198 | 0.956 | 1 | 1 | 2.910 | 1000 |
| pca | ds.vertPCA | eigenvalue_PC2 | 1000 | 0.00294518971 | 0.05878306 | TRUE | 1.000 | 1 | 3.9988830 | 0.956 | 1 | 1 | 2.910 | 1000 |
| pca | ds.vertPCA | eigenvalue_PC3 | 1000 | 0.04956171208 | 0.09017555 | TRUE | 1.000 | 1 | 3.9931663 | 0.956 | 1 | 1 | 2.910 | 1000 |
| pca | ds.vertPCA | PC1_sine_angle | 1000 | 0.06375313469 | 0.08177183 | FALSE | NA | NA | NA | 0.956 | 1 | 1 | 2.910 | 1000 |
| survival | ds.vertDPRMST | rmst_tau_6_fixed_grid | 1000 | 0.00198384050 | 0.05892305 | TRUE | 1.000 | 1 | 1.8685539 | 0.970 | 1 | 1 | 5.221 | 1000 |
| survival | ds.vertDPSurvival | cumulative_incidence_fixed_grid | 12000 | 0.00076107998 | 0.01267566 | TRUE | 1.000 | 1 | 0.2676045 | 0.970 | 1 | 1 | 5.221 | 1000 |
| survival | ds.vertDPSurvival | kaplan_meier_fixed_grid | 6000 | -0.00152215995 | 0.01633827 | TRUE | 1.000 | 1 | 0.4129774 | 0.970 | 1 | 1 | 5.221 | 1000 |
| survival | ds.vertDPSurvival | nelson_aalen_fixed_grid | 6000 | 0.00207497350 | 0.01993186 | TRUE | 1.000 | 1 | 0.5432342 | 0.970 | 1 | 1 | 5.221 | 1000 |
| survival | ds.vertDPSurvivalQuantile | event_probability_0.01_first_grid_endpoint | 1000 | 0.00000000000 | 0.00000000 | TRUE | 1.000 | 1 | 0.0430000 | 0.970 | 1 | 1 | 5.221 | 1000 |

`mechanism_coverage` is conditional on the fixed finite dataset and
covers only DP mechanism noise plus the explicitly reported public-grid
or quantization component. It is not a sampling-coverage estimate.
Gaussian coefficients and the PC1 loading-angle diagnostic are marked
as point-only because dsVert does not claim nonlinear coefficient or
loading regions; their coverage and width fields are therefore `NA`.
The minimum recorded fixed-seed mechanism coverage is `0.994`; this Monte Carlo diagnostic is reported without replacing the analytic and deterministic coverage gates above.

## Boundary cases

| case | pass | detail |
| --- | --- | --- |
| describe_empty_histogram | TRUE | fixed_public_clamp_applied_preclamp_state_not_released |
| quantile_probabilities_0_and_1 | TRUE | dp_projected_histogram_empty,dp_projected_histogram_empty,dp_projected_histogram_empty |
| quantile_nonempty_endpoint_probabilities | TRUE | bin_index=1,8 |
| quantile_exact_cdf_tie_first_crossing | TRUE | bin_index=2,4,6 |
| median_empty_histogram | TRUE | dp_projected_histogram_empty |
| epi_empty_groups | TRUE | dp_point_non_estimable |
| diagnostic_empty_table | TRUE | non_estimable |
| direct_empty_strata | TRUE | dp_point_non_estimable |
| causal_empty_stratum_arms | TRUE | dp_point_non_estimable |
| indirect_zero_expected_denominator | TRUE | dp_point_non_estimable |
| diagnostic_perfect_boundary | TRUE | boundary |
| indirect_infinite_boundary | TRUE | boundary_infinite |
| survival_empty_curve | TRUE | dp_curve_empty_after_postprocessing |
| rmst_empty_curve | TRUE | rmst=6 |
| survival_quantile_beyond_grid | TRUE | mechanism_region_extends_beyond_grid |
| median_survival_beyond_grid | TRUE | mechanism_region_extends_beyond_grid |

## Interpretation

The run exercises the real dsVertClient postprocessors and the deployed
signed-vector radius calculation. Noise draws come from a base-R sampler
for the same ideal two-sided-geometric convolution distribution; they do
not exercise HMAC/HKDF/ChaCha20, Ring128, peer pinning, DSI transport,
sticky replay, server admission, or signed receipt verification. Those
remain separate protocol/E2E validation obligations.
