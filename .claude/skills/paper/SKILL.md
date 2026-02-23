---
name: paper
description: Enter paper revision mode for the SuSiNE LaTeX manuscript. Use when working on writing, figures, or narrative structure.
user-invocable: true
---

# Paper Revision Mode

You are helping revise the SuSiNE paper for Bioinformatics.

## First Steps

Read these files to load context:
1. `refs/SuSiNE - rough draft.tex` (current draft)
2. `refs/susine_study_plan_high_level.md` (what the paper should claim)
3. `refs/susie_susine_background.md` (technical background)

## Draft vs Study Plan Gap

### What exists in the current draft (dated 12/19/25):
- Abstract with 12,600-dataset claim, M1 metric, CBX8 case study
- Model specification (equations 1-4 are correct, keep these)
- Simulation framework section (partially outdated)
- Illustrative example on annotation quality dilution
- CBX8 real-data case study with clustering analysis

### What the study plan requires but is NOT yet in the draft:
- **Ensemble/aggregation framework** (Phase B): exploration mechanisms (random restarts vs prior variance grids vs c-grid) + aggregation mechanisms (ELBO-softmax, uniform, max-ELBO, diagonal BMA)
- **Decomposition of WHY SuSiNE helps** (Phase C): diversity, best-mode quality, weighting/recognition — three axes
- **Failure-mode characterization** (cross-cutting pillar): two-axis map using M1 (LD-only) x Z-primary (z-score-only), candidate z-metric family screening, actionable pre-fit guidance
- **Compute-confound pilot** (Section 5 of plan): random restarts vs c-grid at matched budget K, diversity metrics comparison
- **Multi-locus real data** (Phase D): 3-4 loci across difficulty spectrum (only CBX8 done so far)
- **Explicit hypothesis testing**: H1-H6 need structured reporting

## Hypotheses H1-H6

- **H1**: Multimodality and instability are practical, not pathological edge cases
- **H2**: Ensemble workflows improve SuSiE in a predictable way
- **H3**: SuSiNE's c-grid provides structured exploration that random restarts cannot replicate (D_c-grid >> D_restart at matched compute)
- **H4**: Gains are not "just more compute" — c-grid beats restarts at equal budget K
- **H5**: Directional priors are most valuable when fine-mapping is hardest (high LD complexity, diffuse signal)
- **H6**: Failure modes are predictable from LD-only + z-score-only diagnostics (M1 + z-metric explain operational failure endpoints)

## Writing Guidelines

- **Format**: Bioinformatics style, ~7 pages, ~4 main figures
- **Citations**: Use `\citet{}` (textual) and `\citep{}` (parenthetical)
- **Notation** (must be consistent throughout):
  - `\alpha_\ell` = per-effect PIPs (L vectors, one per single effect)
  - `\mu_0 = c \cdot a` = prior means (annotation scale * annotation vector)
  - `\sigma_0^2` = prior variance on effect sizes
  - `M_1` = mid-LD energy metric (LD-only difficulty diagnostic)
  - JSD = Jensen-Shannon divergence (posterior diversity metric)
  - `p^*` = number of true causal SNPs
  - `\phi_a` = annotation quality (R^2 among causal SNPs)
  - `\phi_y` = noise fraction
- **Methods**: Keep concise; derivations go in supplement
- **"Use case" terminology**: Maps to model configurations from use_case_catalog() in R/use_cases.R (14 configs: a_i through c_ii)
- **Tone**: Precise, quantitative, suitable for computational biology audience