# Catalogs -------------------------------------------------------------------

#' Prior-spec catalog (model definitions only).
#'
#' Each row is a model specification independent of exploration/aggregation.
#'
#' @return Tibble with prior-spec metadata.
#' @export
prior_spec_catalog <- function() {
  tibble::tribble(
    ~prior_spec_id,                         ~use_case_id,                          ~label,                                                    ~backend,  ~prior_mean_strategy, ~prior_variance_strategy, ~inclusion_prior_strategy, ~unmappable_effects, ~eb_method,                ~c_nonneg,
    # ---- Existing pilot use cases -------------------------------------------------------
    "susie_vanilla",                        "susie_vanilla",                        "SuSiE vanilla (susieR)",                                  "susieR",  "zero",               "fixed",                  "uniform",                 "none",              NA_character_,             FALSE,
    "susie_ash_fixed",                      "susie_ash_fixed",                      "SuSiE-ash fixed prior var (susieR)",                      "susieR",  "zero",               "fixed",                  "uniform",                 "ash",               NA_character_,             FALSE,
    "susie_ash_eb",                         "susie_ash_eb",                         "SuSiE-ash EB prior var (susieR)",                         "susieR",  "zero",               "eb",                     "uniform",                 "ash",               NA_character_,             FALSE,
    "susie_inf",                            "susie_inf",                            "SuSiE-inf (susieR)",                                      "susieR",  "zero",               "fixed",                  "uniform",                 "inf",               NA_character_,             FALSE,
    "susine_vanilla",                       "susine_vanilla",                       "SuSiE vanilla (susine)",                                  "susine",  "zero",               "fixed",                  "uniform",                 "none",              NA_character_,             FALSE,
    "susie_eb",                             "susie_eb",                             "SuSiE EB prior var (susieR)",                             "susieR",  "zero",               "eb",                     "uniform",                 "none",              NA_character_,             FALSE,
    "susine_eb",                            "susine_eb",                            "SuSiE EB prior var (susine)",                             "susine",  "zero",               "eb",                     "uniform",                 "none",              "var",                     FALSE,
    "susie_functional_pi",                  "susie_functional_pi",                  "SuSiE functional pi (susieR)",                            "susieR",  "zero",               "fixed",                  "functional_pi",           "none",              NA_character_,             FALSE,
    "susine_functional_pi",                 "susine_functional_pi",                 "SuSiE functional pi (susine)",                            "susine",  "zero",               "fixed",                  "functional_pi",           "none",              NA_character_,             FALSE,
    "susine_functional_mu",                 "susine_functional_mu",                 "SuSiNE functional mu_0 = c*a (susine)",                   "susine",  "functional_mu",      "fixed",                  "uniform",                 "none",              NA_character_,             FALSE,
    # ---- EB method ablation use cases (all: susine backend, functional_mu, uniform pi) --
    # var: EB on sigma only; annotation is set as mu_0=c*a but c stays fixed at c_value
    "susine_eb_annotated",                  "susine_eb_annotated",                  "SuSiNE EB var (annotation, susine)",                      "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "var",                     FALSE,
    # mean: EB on scalar mean only; sigma fixed; no annotation needed
    "susine_eb_mean",                       "susine_eb_mean",                       "SuSiNE EB mean only (susine)",                            "susine",  "zero",               "eb",                     "uniform",                 "none",              "mean",                    FALSE,
    # mu_var: EB on scalar mean + sigma jointly; no annotation needed
    "susine_eb_mu_var",                     "susine_eb_mu_var",                     "SuSiNE EB mean+var (susine)",                             "susine",  "zero",               "eb",                     "uniform",                 "none",              "mu_var",                  FALSE,
    # scale: EB on c only (sigma fixed); c_nonneg=FALSE
    "susine_eb_scale",                      "susine_eb_scale",                      "SuSiNE EB c only (susine)",                               "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "scale",                   FALSE,
    # scale: EB on c only (sigma fixed); c_nonneg=TRUE
    "susine_eb_scale_nonneg",               "susine_eb_scale_nonneg",               "SuSiNE EB c only nonneg (susine)",                        "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "scale",                   TRUE,
    # scale_var: joint EB on c + sigma; c_nonneg=FALSE
    "susine_eb_scale_var",                  "susine_eb_scale_var",                  "SuSiNE EB c+var (susine)",                                "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "scale_var",               FALSE,
    # scale_var: joint EB on c + sigma; c_nonneg=TRUE
    "susine_eb_scale_var_nonneg",           "susine_eb_scale_var_nonneg",           "SuSiNE EB c+var nonneg (susine)",                         "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "scale_var",               TRUE,
    # penalized_scale_var: joint EB with IG penalty on sigma; c_nonneg=FALSE
    "susine_eb_penalized_scale_var",        "susine_eb_penalized_scale_var",        "SuSiNE EB penalized c+var (susine)",                      "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "penalized_scale_var",     FALSE,
    # penalized_scale_var: joint EB with IG penalty on sigma; c_nonneg=TRUE
    "susine_eb_penalized_scale_var_nonneg", "susine_eb_penalized_scale_var_nonneg", "SuSiNE EB penalized c+var nonneg (susine)",               "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "penalized_scale_var",     TRUE,
    # clamped_scale_var: joint EB with sigma floor; c_nonneg=FALSE
    "susine_eb_clamped_scale_var",          "susine_eb_clamped_scale_var",          "SuSiNE EB clamped c+var (susine)",                        "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "clamped_scale_var",       FALSE,
    # clamped_scale_var: joint EB with sigma floor; c_nonneg=TRUE
    "susine_eb_clamped_scale_var_nonneg",   "susine_eb_clamped_scale_var_nonneg",   "SuSiNE EB clamped c+var nonneg (susine)",                 "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "clamped_scale_var",       TRUE,
    # alternating_scale_var: alternating 1D Brent steps; c_nonneg=FALSE
    "susine_eb_alternating_scale_var",      "susine_eb_alternating_scale_var",      "SuSiNE EB alternating c+var (susine)",                    "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "alternating_scale_var",   FALSE,
    # alternating_scale_var: alternating 1D Brent steps; c_nonneg=TRUE
    "susine_eb_alternating_scale_var_nonneg","susine_eb_alternating_scale_var_nonneg","SuSiNE EB alternating c+var nonneg (susine)",            "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "alternating_scale_var",   TRUE,
    # scale_var_outer: inner IBSS does var only, outer loop re-estimates c; c_nonneg=FALSE
    "susine_eb_scale_var_outer",            "susine_eb_scale_var_outer",            "SuSiNE EB outer c+var (susine)",                          "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "scale_var_outer",         FALSE,
    # scale_var_outer: inner IBSS does var only, outer loop re-estimates c; c_nonneg=TRUE
    "susine_eb_scale_var_outer_nonneg",     "susine_eb_scale_var_outer_nonneg",     "SuSiNE EB outer c+var nonneg (susine)",                   "susine",  "functional_mu",      "eb",                     "uniform",                 "none",              "scale_var_outer",         TRUE
  ) %>%
    dplyr::mutate(
      supports_annotation = .data$prior_mean_strategy == "functional_mu" |
        .data$inclusion_prior_strategy == "functional_pi"
    )
}

#' Exploration catalog.
#'
#' @return Tibble of exploration methods.
#' @export
exploration_catalog <- function() {
  tibble::tribble(
    ~exploration_id, ~label,
    "single",  "Single fit",
    "restart", "Restart via randomized prior inclusion vectors",
    "c_grid",  "Grid over c for mu_0 = c*a",
    "tau_grid","Grid over tau for pi = softmax(|a|/tau)",
    "refine",        "Harness-level refinement (perturb + refit BFS)",
    "sigma_0_2_grid","Grid over prior variance sigma_0_2"
  )
}

#' Aggregation catalog.
#'
#' @return Tibble of aggregation methods.
#' @export
aggregation_catalog <- function() {
  tibble::tribble(
    ~aggregation_id, ~label,
    "max_elbo",          "Max ELBO",
    "uniform",           "Uniform average",
    "elbo_softmax",      "ELBO softmax",
    "cluster_weight",        "Cluster-then-ELBO-softmax (JSD 0.15)",
    "cluster_weight_jsd_050","Cluster-then-ELBO-softmax (JSD 0.50)"
  )
}

#' Validate exploration compatibility with prior specs.
#'
#' @param prior_spec_ids Character vector of prior_spec_id values.
#' @param exploration_ids Character vector of exploration_id values.
#' @return Tibble with columns prior_spec_id, exploration_id, valid.
#' @export
valid_exploration_for_prior <- function(prior_spec_ids = NULL,
                                        exploration_ids = NULL) {
  prior_specs <- prior_spec_catalog()
  explorations <- exploration_catalog()
  if (!is.null(prior_spec_ids)) {
    prior_specs <- dplyr::filter(prior_specs, .data$prior_spec_id %in% prior_spec_ids)
  }
  if (!is.null(exploration_ids)) {
    explorations <- dplyr::filter(explorations, .data$exploration_id %in% exploration_ids)
  }
  tidyr::expand_grid(
    prior_spec_id = prior_specs$prior_spec_id,
    exploration_id = explorations$exploration_id
  ) %>%
    dplyr::left_join(
      dplyr::select(
        prior_specs,
        .data$prior_spec_id,
        .data$prior_mean_strategy,
        .data$prior_variance_strategy,
        .data$inclusion_prior_strategy,
        .data$eb_method
      ),
      by = "prior_spec_id"
    ) %>%
    dplyr::mutate(
      # EB methods that optimize c (annotation scale) — c_grid is redundant
      .eb_optimizes_c = !is.na(.data$eb_method) & .data$eb_method %in% c(
        "scale", "scale_var", "scale_var_outer",
        "penalized_scale_var", "clamped_scale_var", "alternating_scale_var"
      ),
      valid = dplyr::case_when(
        .data$exploration_id == "single" ~ TRUE,
        .data$exploration_id == "restart" ~ TRUE,
        .data$exploration_id == "refine" ~ TRUE,
        .data$exploration_id == "c_grid" ~
          .data$prior_mean_strategy == "functional_mu" & !.data$.eb_optimizes_c,
        .data$exploration_id == "tau_grid" ~ .data$inclusion_prior_strategy == "functional_pi",
        .data$exploration_id == "sigma_0_2_grid" ~
          .data$prior_variance_strategy == "fixed",
        TRUE ~ FALSE
      )
    ) %>%
    dplyr::select(-.data$.eb_optimizes_c) %>%
    dplyr::select(.data$prior_spec_id, .data$exploration_id, .data$valid)
}

#' Backward-compatible use-case catalog alias.
#'
#' Returns the prior-spec catalog with `use_case_id` naming.
#'
#' @return Tibble model catalog.
#' @export
use_case_catalog <- function() {
  prior_spec_catalog() %>%
    dplyr::select(
      .data$use_case_id,
      .data$prior_spec_id,
      .data$label,
      .data$backend,
      .data$prior_mean_strategy,
      .data$prior_variance_strategy,
      .data$inclusion_prior_strategy,
      .data$unmappable_effects,
      .data$supports_annotation
    )
}

#' Resolve model specs by id.
#'
#' Accepts either prior_spec_id values or use_case_id aliases.
#'
#' @param ids Character vector.
#' @return Tibble subset of prior specs.
#' @export
resolve_use_cases <- function(ids = NULL) {
  catalog <- prior_spec_catalog()
  if (is.null(ids)) {
    return(catalog)
  }
  dplyr::filter(catalog, .data$prior_spec_id %in% ids | .data$use_case_id %in% ids)
}
