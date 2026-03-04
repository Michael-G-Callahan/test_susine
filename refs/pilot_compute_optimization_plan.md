# Pilot Study Compute Optimization Plan

*Created 2026-03-04. Context: pilot_study_v1 overran budget.*

## Context

The pilot study (`pilot_study_v1`) ran 90 SLURM array tasks for ~15 hours, burning ~2.79 credits (vs. 0.0496 estimated). Only 20-30 of 90 tasks produced even one flush. The job has **1,648,350 total runs** -- far too many. Goal: complete a meaningful pilot for <1 credit total, using a two-phase approach: screen use cases cheaply first, then explore only the promising ones.

---

## Diagnostics

### Why ~2.79 credits instead of 0.0496?

`job_estimate` computed cost for **one array element**, not the full array. Multiply by 90: 0.0496 x 90 = 4.46 credits (upper bound at full 24h). Actual burn was ~2.79 because tasks were cancelled before completing.

### Wall-time patterns

Per-run times are mostly fine (median 0.68s, mean 1.3s). The problem is volume (1.65M runs). Outliers: `refine` has a 60s worst case, `restart` has 46s worst cases -- but rare (<1% of runs).

### Run accounting -- where the compute goes

**91% of runs are in 3 functional use cases** (502K each: susie_functional_pi, susine_functional_pi, susine_functional_mu). Each crosses: 48 groups/bundle x (1 single + 10 restarts + 10 refine + grid values) x 3 sigma_0_2 x 450 bundles.

### Memory/flush

All ~3,663 runs per bundle accumulate in memory before flushing. With the lean design below (~30-150 runs/bundle), this concern disappears.

---

## Two-Phase Pilot Design

### Phase A: Screening (No Exploration)

**Purpose**: All 10 use cases, single fits only. Screen which use cases merit exploration.

**Grid**:

- **Use cases**: All 10
- **Exploration**: `single` only
- **Matrices**: 10 (M1-stratified)
- **p_star**: 5 levels (1, 3, 5, 10, 20)
- **y_noise**: 3 levels (0.5, 0.8, 0.95)
- **annotation_r2**: 3 levels (0.0, 0.1, 0.5)
- **inflate_match**: 3 levels (0.8, 1.0, 1.2)
- **sigma_0_2**: 3 levels (0.1, 0.2, 0.4) for fixed; 1 for EB
- **Seeds**: 3

**Estimated runs**: ~46,000

**SLURM**: `--mem=6G`, `--time=23:59:59`, `bundles_per_task = 10` (~45 tasks)

**Decision criteria**:

- If EB use cases underperform baselines -> drop from Phase B
- If susie_functional_pi ~ susine_functional_pi -> keep only one
- Identify uninformative grid dimensions for trimming

### Phase B: Exploration (Targeted)

**Purpose**: Add exploration to Phase-A survivors.

**Expected use cases** (adjust based on Phase A):

1. `susie_vanilla` -- baseline anchor
2. `susie_inf` or `susie_ash_fixed` -- best non-SuSiNE competitor
3. `susine_functional_mu` -- central method (c-grid)
4. `susine_functional_pi` -- functional pi arm (tau-grid)

**Exploration**:

- Baselines: `single` + `restart` (K=5)
- `susine_functional_mu`: `single` + `c_grid` (K=5) + `restart` (K=5)
- `susine_functional_pi`: `single` + `tau_grid` (K=5) + `restart` (K=5)

**Grid**: Same as Phase A, potentially trimmed. sigma_0_2: 2 levels for functional (0.1, 0.3); 1 (0.2) for baselines.

**Estimated runs**: ~280,000

**Combined budget**: Well under 1 credit.

---

## Workbooks

- Phase A: `vignettes/simulation pipeline/run_control_workbook_phase_a.Rmd`
- Phase B: `vignettes/simulation pipeline/run_control_workbook_phase_b.Rmd`
- Original: `vignettes/simulation pipeline/run_control_workbook_pilot.Rmd`
