# Catalogs -------------------------------------------------------------------

#' Prior-spec catalog (model definitions only).
#'
#' Each row is a model specification independent of exploration/aggregation.
#'
#' @return Tibble with prior-spec metadata.
#' @export
prior_spec_catalog <- function() {
  tibble::tribble(
    ~prior_spec_id, ~use_case_id, ~label, ~backend, ~prior_mean_strategy, ~prior_variance_strategy, ~inclusion_prior_strategy, ~unmappable_effects,
    "susie_vanilla",      "susie_vanilla",       "SuSiE vanilla (susieR)",               "susieR", "zero",          "fixed", "uniform",      "none",
    "susie_ash",          "susie_ash",           "SuSiE-ash (susieR)",                   "susieR", "zero",          "fixed", "uniform",      "ash",
    "susie_inf",          "susie_inf",           "SuSiE-inf (susieR)",                   "susieR", "zero",          "fixed", "uniform",      "inf",
    "susine_vanilla",     "susine_vanilla",      "SuSiE vanilla (susine)",               "susine", "zero",          "fixed", "uniform",      "none",
    "susie_eb",           "susie_eb",            "SuSiE EB prior var (susieR)",          "susieR", "zero",          "eb",    "uniform",      "none",
    "susine_eb",          "susine_eb",           "SuSiE EB prior var (susine)",          "susine", "zero",          "eb",    "uniform",      "none",
    "susie_functional_pi","susie_functional_pi", "SuSiE functional pi (susieR)",         "susieR", "zero",          "fixed", "functional_pi","none",
    "susine_functional_pi","susine_functional_pi","SuSiE functional pi (susine)",        "susine", "zero",          "fixed", "functional_pi","none",
    "susine_functional_mu","susine_functional_mu","SuSiNE functional mu_0 = c*a (susine)","susine","functional_mu", "fixed", "uniform",      "none"
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
    "refine",  "Harness-level refinement (BFS)"
  )
}

#' Aggregation catalog.
#'
#' @return Tibble of aggregation methods.
#' @export
aggregation_catalog <- function() {
  tibble::tribble(
    ~aggregation_id, ~label,
    "max_elbo",      "Max ELBO",
    "uniform",       "Uniform average",
    "elbo_softmax",  "ELBO softmax",
    "cluster_weight","Cluster-then-ELBO-softmax (Method A)"
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
        .data$inclusion_prior_strategy
      ),
      by = "prior_spec_id"
    ) %>%
    dplyr::mutate(
      valid = dplyr::case_when(
        .data$exploration_id == "single" ~ TRUE,
        .data$exploration_id == "restart" ~ TRUE,
        .data$exploration_id == "refine" ~ TRUE,
        .data$exploration_id == "c_grid" ~ .data$prior_mean_strategy == "functional_mu",
        .data$exploration_id == "tau_grid" ~ .data$inclusion_prior_strategy == "functional_pi",
        TRUE ~ FALSE
      )
    ) %>%
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
