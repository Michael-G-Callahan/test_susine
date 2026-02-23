---
name: review-r
description: Review R code for correctness, style, and performance. Tailored to this statistical research codebase.
user-invocable: true
---

# R Code Review

Review the provided R code using this project-specific checklist. Report findings organized by severity: Bugs > Risks > Style > Performance.

## Correctness

- [ ] **NA handling**: Are NAs propagated safely? Check `na.rm` usage in `mean()`, `sum()`, `max()`. Watch for silent NA drops in joins.
- [ ] **Numerical stability**: Check for `log(0)`, division by zero, underflow in ELBO/BF/softmax computations. Use `log1p()` and `expm1()` where appropriate.
- [ ] **Edge cases**: Zero-length vectors, single-row tibbles, all-NA columns, L=1, p_star=0.
- [ ] **Seed management**: Is `set.seed()` used correctly? Could parallel execution cause seed collisions?
- [ ] **Type safety**: Watch for `as.integer(NA)` vs `NA_integer_`, `TRUE/FALSE` vs `1L/0L`, factor vs character coercion.
- [ ] **Matrix dimensions**: Verify dimensions match before matrix multiplication. Check drop=FALSE on single-column subsets.
- [ ] **Vectorization correctness**: Ensure vectorized operations broadcast as intended (recycling rules).

## Style (project conventions)

- [ ] Assignment: `<-`, not `=`
- [ ] Pipe: `%>%`, not `|>`
- [ ] Tidy eval: `.data$col` in dplyr verbs
- [ ] NULL coalescing: `%||%` for defaults (defined in utils.R)
- [ ] Naming: `snake_case` for all identifiers
- [ ] roxygen: Markdown format, `@keywords internal` for non-exported helpers
- [ ] Section headers: `# Section name -----`
- [ ] Tibbles over data.frames
- [ ] Explicit namespace: `dplyr::filter()` not bare `filter()` in package code

## Performance (relevant for HPC at scale: 12,600+ datasets x 77 configs)

- [ ] **Memory**: Are large objects (matrices, parquet) freed with `rm()` + `gc()` after use?
- [ ] **Allocation**: Pre-allocate vectors in loops instead of growing via `c()` or `rbind()`.
- [ ] **I/O**: Use `arrow::write_parquet()` for large tables, not CSV.
- [ ] **Caching**: Could repeated computations (e.g., `R = cor(X)`, `standardize_x()`) be cached?
- [ ] **Matrix ops**: Use `crossprod(X)` instead of `t(X) %*% X`. Use `tcrossprod()` where appropriate.
- [ ] **Column stats**: Prefer `matrixStats::colMeans2()` / `colVars()` over `apply(X, 2, ...)`.

## Domain-Specific (statistical correctness)

- [ ] **PIP computation**: Formula should be `1 - prod(1 - alpha_l)` across L effects.
- [ ] **Credible sets**: Sorted by decreasing alpha, cumulative sum >= rho threshold.
- [ ] **CS purity**: Min absolute pairwise correlation in the set, not mean.
- [ ] **AUPRC**: Trapezoidal integration, anchored at (recall=0, precision=1).
- [ ] **ELBO**: Should be monotonically non-decreasing across iterations (up to numerical noise).
- [ ] **Softmax weighting**: `exp(ELBO - max(ELBO))` for numerical stability, then normalize.
- [ ] **JSD**: Symmetric, bounded [0, 1] when using log base 2. Check `M = 0.5*(P+Q)`.

## Output Format

```
## Bugs (incorrect logic)
- [file:line] Description of the bug and its impact

## Risks (fragile code)
- [file:line] What could go wrong and under what conditions

## Style (convention violations)
- [file:line] What convention is violated

## Performance (optimization opportunities)
- [file:line] What could be faster and estimated impact
```