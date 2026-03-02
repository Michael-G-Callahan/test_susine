# test_susine Package Updates: Detailed Scope and Refactor Plan

**Date:** 2026-02-27  
**Status:** Planning only (no implementation in this document)

## 1) Scope and Goal

This plan scopes all `T1`-`T20` items from `refs/project_status.md` Section 3.2, maps each item to the current codebase, and orders implementation for lowest integration risk.

Primary objective:
- Move from the legacy use-case-centric harness (single `susine` backend + bundled exploration logic) to a factored design:
  - prior spec catalog
  - exploration catalog
  - aggregation catalog
  - explicit compatibility constraints

## 2) Current-State Findings (from code)

1. The model runner currently dispatches only to `susine::susine`.
   - Evidence: `R/run_model.R` calls `susine::susine` in `run_use_case()`.
2. Use cases are still the legacy 14-row catalog (`a_i ... c_ii`) with exploration mixed into use-case rows (`extra_compute`).
   - Evidence: `R/use_cases.R`.
3. Run-control expansion is use-case-driven, not factored by independent prior/exploration/aggregation catalogs.
   - Evidence: `R/run_controls.R` (`make_run_tables()`, `make_job_config()`).
4. Aggregation methods currently include `mean`, `max_elbo`, `softmax_elbo`; Method A cluster-then-weight is not implemented.
   - Evidence: `R/run_model.R` (`aggregate_use_case_pips()`).
5. Multimodal clustering currently uses average linkage; Method A spec recommends complete linkage.
   - Evidence: `R/run_model.R` (`compute_multimodal_metrics()`).
6. Simulation architecture is sparse-only; no oligogenic architecture flag/path exists.
   - Evidence: `R/simulate_data.R` (`simulate_effect_sizes()`).
7. `gamma_shrink` is still active in simulation and run-control grids (not deprecated).
   - Evidence: `R/simulate_data.R`, `R/run_controls.R`, `R/run_model.R`.
8. Per-fit PIP vectors are already persisted in parquet (`snps.parquet` / flush shards), and ELBO final is in model metrics.
   - Evidence: `R/run_model.R` (`build_snp_table()`, `write_run_outputs()`).
9. Wall-time per fit is not currently recorded.
   - Evidence: no timing instrumentation around fit calls in `R/run_model.R`.
10. Dataset metrics include `M1` + z metrics, but not explicit off-diagonal `|r| > 0.95` count metric.
    - Evidence: `R/dataset_metrics.R` (`compute_dataset_metrics()`).
11. Seed handling needs audit: `dataset_bundles$phenotype_seed` is overwritten by row number after bundle construction.
    - Evidence: `R/run_controls.R` (`make_run_tables()`).

## 3) Coupling and Dependency Rules

Hard coupling:
- `T1 + T2 + T12` must be implemented as one coherent refactor unit.

Strong dependencies:
- `T5`, `T8`, `T11`, `T16`, `T17` depend on the refactored catalogs/control plane from `T1+T2+T12`.
- `T6` depends on per-fit persistence and aggregation hooks from the same refactor.
- `T19` depends on `T6` (or can be finalized during `T6`).
- `T10` depends on factored aggregation controls and composable exploration outputs.
- `T20` should run after architecture+wiring are stable.

Cross-cutting validation dependencies:
- `T13` (metrics audit) should run after the first integrated pass of `T1+T2+T12+T5+T6+T16+T17`.
- `T18` (seed management review) should be pulled earlier than currently listed because it affects reproducibility of every subsequent test run.

## 4) Recommended Refactor Order

## Phase A: Architecture Refactor (single branch block)
1. `T1 + T2 + T12` (single coupled change)
2. `T18` seed management review (early, before broad testing)

## Phase B: Wiring and Core Feature Completion
3. `T8` expose exploration budget `K` as top-level workbook/control parameter
4. `T5` tau-grid exploration for functional `pi`
5. `T11` exploration composability (e.g., restart x c-grid x tau-grid)
6. `T16` softmax temperature threading from job config to aggregation
7. `T17` per-fit wall-time recording

## Phase C: Aggregation and Metrics
8. `T6` Method A cluster-then-ELBO-softmax
9. `T19` finalize complete-vs-average linkage decision
10. `T10` overall per-dataset aggregation across all fits/use-cases
11. `T13` full metrics persistence audit
12. `T14` high-LD-count metric
13. `T15` verify/confirm per-fit PIP persistence contract (mostly already present)

## Phase D: Simulation and Legacy Cleanup
14. `T3` oligogenic architecture
15. `T9` deprecate `gamma_shrink`
16. `T7` default vs warm start labeling
17. `T4` refinement exploration mode (BFS, model-agnostic)

## Phase E: Integration Validation
18. `T20` 2-locus end-to-end smoke test

Rationale for this order:
- Architecture first avoids repeated rewrites.
- Seed control is moved earlier to prevent invalid/irreproducible validation.
- Exploration/aggregation wiring precedes advanced simulation additions.
- BFS refinement is intentionally delayed until the new control plane is stable.

## 5) Detailed Plan by Task

## T1 - Wire susieR use cases into harness
Scope:
- Add backend dispatch layer in `R/run_model.R` so runner can call either:
  - `susine::susine`
  - `susieR::susie`
- Introduce a backend-specific argument mapper:
  - `prior_inclusion_weights` -> `prior_weights`
  - `sigma_0_2` -> `scaled_prior_variance`
  - `init_alpha`/warm starts -> `model_init` strategy
  - `eb_prior_mode` equivalents -> `estimate_prior_variance`
- Normalize fit outputs into a backend-agnostic internal structure used by `evaluate_model()` and write paths.

Touched files:
- `R/run_model.R`
- `R/use_cases.R` (or replacement catalogs)
- `DESCRIPTION`/`NAMESPACE` only if hard import strategy changes

Done criteria:
- Same run-control can execute both susieR and susine prior specs.
- Output schema remains stable for downstream aggregation/collection.

## T2 - Refactor use case catalog
Scope:
- Replace monolithic `use_case_catalog()` with:
  - `prior_spec_catalog()`
  - `exploration_catalog()`
  - `aggregation_catalog()`
  - `valid_exploration_for_prior()`
- Remove old 14-entry catalog without compatibility shim.

Touched files:
- `R/use_cases.R` (replace with catalog module)
- `R/run_controls.R`
- `R/visualize_results.R` (label lookup adjustments)
- vignettes in `vignettes/simulation pipeline/` (input schema changes)

Done criteria:
- Catalogs are independent and fully drive run construction.
- All selected prior specs from Section 2.1 are representable.

## T3 - Implement oligogenic architecture
Scope:
- Add `architecture` control in simulation spec.
- Implement `simulate_effect_sizes_oligogenic()` with tiered h2 partitioning.
- Ensure annotation simulation path is architecture-aware where required.

Touched files:
- `R/simulate_data.R`
- `R/run_controls.R` (new grid knobs)
- run-control workbooks

Done criteria:
- Sparse and oligogenic simulations can be selected through run config.
- Metadata clearly records architecture in outputs.

## T4 - Implement refinement exploration mode (BFS)
Scope:
- Add model-agnostic harness-level BFS exploration.
- Do not rely on susieR `refine=TRUE` for this task.
- Preserve all generated fits for ensembling (no greedy discard).

Touched files:
- `R/run_model.R`
- `R/run_controls.R` / exploration catalog

Done criteria:
- Exploration method selectable with budget cap `K`.
- Produces fit set compatible with existing aggregation and metrics pipeline.

## T5 - Implement tau-grid exploration for functional pi
Scope:
- Add `tau_grid` to config and expansion logic.
- Wire tau into functional-pi prior spec for both susieR and susine-equivalent specs where applicable.

Touched files:
- `R/run_controls.R`
- `R/run_model.R`
- run-control workbook(s)

Done criteria:
- Functional-pi runs expand over tau values analogously to current c-grid behavior.

## T6 - Implement cluster-then-ELBO-softmax (Method A)
Scope:
- Add Method A aggregation:
  - cluster fits by JSD
  - choose representative per cluster (max ELBO)
  - apply importance-corrected weights `exp(ELBO)/f_m`
  - report ESS

Touched files:
- `R/run_model.R` (aggregation function family)
- possibly `R/collect_results.R` if extra outputs/diagnostics are persisted

Done criteria:
- `aggregation_methods` can include cluster-weight method.
- Outputs include method identifier and diagnostics (at least ESS).

## T7 - Label default vs warm start fits
Scope:
- Add `init_type` in fit metadata: `default` or `warm`.
- For restart groups ensure first fit is explicit default (no warm start init object), remaining are warm.

Touched files:
- `R/run_model.R`
- output schemas in write paths

Done criteria:
- Downstream outputs can stratify results by initialization type.

## T8 - Pass warm start K through run-control workbook
Scope:
- Add explicit top-level `K` exploration budget in workbook and `make_job_config()`.
- Map `K` to applicable exploration methods (restarts, refinement breadth, grid truncation policy as needed).

Touched files:
- `R/run_controls.R`
- `vignettes/simulation pipeline/run_control_workbook*.Rmd`

Done criteria:
- User can set one budget knob in workbook and see it reflected in run table expansion.

## T9 - Deprecate gamma_shrink in simulate_priors()
Scope:
- Add deprecation warning path for `gamma_shrink`.
- Stop emitting `gamma_shrink` in new configs.
- Keep backward compatibility for reading legacy run tables during transition.

Touched files:
- `R/simulate_data.R`
- `R/run_controls.R`
- `R/run_model.R` (cache key path simplification after transition)

Done criteria:
- New runs no longer depend on gamma-shrink behavior.
- Legacy runs remain readable.

## T10 - Add overall per-dataset aggregation option
Scope:
- Add control to pool all fits for a dataset across use-cases/specs before aggregation.
- Keep existing per-group aggregation path alongside global mode.

Touched files:
- `R/run_model.R`
- `R/run_controls.R`

Done criteria:
- Config can toggle per-use-case aggregation, global aggregation, or both.

## T11 - Allow exploration composability
Scope:
- Ensure expansions can compose multiplicatively where valid:
  - restart x c-grid x tau-grid x (optional refine seeds/branches)
- Drive this from exploration catalog rather than hardcoded branching.

Touched files:
- `R/run_controls.R`
- `R/run_model.R`

Done criteria:
- Run table clearly encodes composite exploration provenance for each fit.

## T12 - Separate exploration/aggregation from use cases
Scope:
- Ensure prior specs define model only.
- Exploration and aggregation are independent control inputs.
- Remove `extra_compute` from prior-spec definition.

Touched files:
- `R/use_cases.R` replacement
- `R/run_controls.R`
- `R/run_model.R`

Done criteria:
- No exploration semantics are stored in prior-spec rows.

## T13 - Verify metrics are collected and saved
Scope:
- Audit outputs against Section 2.7 requirements:
  - model metrics
  - CS metrics
  - dataset metrics
  - multimodal metrics
- Verify persistence paths in both verbose and buffered modes.

Touched files:
- mostly audit in `R/run_model.R`, `R/collect_results.R`, `R/evaluation_metrics.R`, `R/dataset_metrics.R`

Done criteria:
- Checklist document maps each required metric to output file/column.
- Any missing metric has implementation issue logged.

## T14 - Add high-LD-count metric
Scope:
- Add off-diagonal count `sum(|r| > 0.95)` (and optional normalized transform).

Touched files:
- `R/dataset_metrics.R`

Done criteria:
- New metric is present in dataset metrics output and documented.

## T15 - Ensure per-fit PIP vectors are saved
Scope:
- Confirm per-fit PIPs and ELBOs are retained for all fit variants.
- Validate aggregation inputs can be reconstructed from stored artifacts.

Current state:
- Largely implemented via SNP parquet + model metrics.

Touched files:
- likely audit only unless schema gaps are found.

Done criteria:
- Explicit guarantee and test checklist for fit-level recoverability.

## T16 - Wire softmax temperature through config
Scope:
- Add `softmax_temperature` to job config.
- Thread into softmax aggregation weight function.

Touched files:
- `R/run_controls.R`
- `R/run_model.R`
- workbook(s)

Done criteria:
- Different temperatures alter aggregation output deterministically.

## T17 - Add wall time recording per fit
Scope:
- Instrument fit execution timing around each model call.
- Persist per-fit timing in model metrics output.

Touched files:
- `R/run_model.R`

Done criteria:
- Every fit record has wall-time field in model metrics.

## T18 - Seed management review
Scope:
- Audit seed propagation for:
  - phenotype generation
  - annotation draws
  - restart/warm-start initialization
  - exploration branch creation
- Define deterministic seed derivation rules and document them.

Touched files:
- `R/run_controls.R`
- `R/run_model.R`
- workbook defaults/documentation

Done criteria:
- Re-running with same config reproduces same data and fit-level outputs.

## T19 - Complete vs average linkage evaluation
Scope:
- Evaluate linkage choice specifically for Method A cluster definition.
- Keep average linkage for multimodal summary only if justified separately.

Touched files:
- `R/run_model.R` (clustering call sites)

Done criteria:
- Linkage policy is explicit and documented.

## T20 - 2-locus pipeline smoke test
Scope:
- Execute full end-to-end small run after refactor integration.
- Validate artifacts and aggregation tables are produced without manual fixes.

Touched files:
- no core changes required; uses run-control + scripts.

Done criteria:
- Smoke-test checklist passes and failures (if any) are logged before pilot launch.

## 6) File-Level Work Packages (for implementation phase)

Package A: Control-plane refactor
- `R/use_cases.R`
- `R/run_controls.R`
- workbooks in `vignettes/simulation pipeline/`

Package B: Runner/backend + fit metadata
- `R/run_model.R`
- `R/evaluation_metrics.R` (if backend normalization requires)

Package C: Metrics and aggregation
- `R/dataset_metrics.R`
- `R/run_model.R`
- `R/collect_results.R`

Package D: Simulation architecture
- `R/simulate_data.R`
- `R/run_controls.R`

Package E: Validation and smoke tests
- scripts/workbooks + aggregation checks

## 7) Risk Register (pre-implementation)

1. API mismatch risk (susieR vs susine output objects).
2. Run-table explosion risk after composable exploration.
3. Reproducibility risk due to current seed handling.
4. Output-schema drift risk affecting existing collection/visualization notebooks.
5. Compute-budget fairness risk if `K` semantics differ by exploration method.

Mitigation:
- Freeze/define internal fit adapter schema first.
- Add run-count accounting before and after refactor.
- Lock deterministic seeding rules early.
- Keep schema migration notes and version stamps in job config.

