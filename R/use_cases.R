# Use case catalog -----------------------------------------------------------

#' Predefined model configurations for simulation sweeps.
#'
#' Each row corresponds to one scenario described in the project brief.
#'
#' @return A tibble with metadata describing each use case.
#' @export
use_case_catalog <- function() {
  tibble::tribble(
    ~use_case_id, ~group, ~label, ~mu_strategy, ~sigma_strategy,
    ~prior_update_method, ~auto_scale_mu, ~auto_scale_sigma,
    ~extra_compute, ~requires_prior_quality,
    "a_i",   "1a", "SuSiE - naive sigma, naive mu",           "naive",      "naive",
    "none",       FALSE, FALSE, "none",       FALSE,
    "a_ii",  "1a", "SuSiE - EB sigma, naive mu",              "naive",      "eb_sigma",
    "var",        FALSE, FALSE, "none",       FALSE,
    "a_iii", "1a", "SuSiNE - EB mu, naive sigma",             "eb_mu",      "naive",
    "mean",       FALSE, FALSE, "none",       FALSE,
    "a_iv",  "1a", "SuSiNE - EB mu & sigma",                  "eb_mu",      "eb_sigma",
    "both",       FALSE, FALSE, "none",       FALSE,
    "b_i",   "1b", "SuSiE + functional sigma (mu=0)",        "naive",      "functional",
    "none",       FALSE, FALSE, "none",       TRUE,
    "b_ii",  "1b", "SuSiNE + functional mu (sigma naive)",    "functional", "naive",
    "none",       FALSE, FALSE, "none",       TRUE,
    "b_iii", "1b", "SuSiNE + functional mu (sigma EB)",       "functional", "eb_sigma",
    "none",       FALSE, FALSE, "none",       TRUE,
    "c_i",   "1c", "SuSiNE + tempering/annealing",            "naive",      "naive",
    "none",       FALSE, FALSE, "anneal",     FALSE,
    "c_ii",  "1c", "SuSiNE + model averaging (multi-init)",   "naive",      "naive",
    "none",       FALSE, FALSE, "model_avg",  FALSE
  )
}

#' Fetch use-case metadata for a given identifier vector.
#'
#' @param ids Character vector of use_case_id entries.
#' @return Tibble subset of `use_case_catalog()`.
#' @export
resolve_use_cases <- function(ids) {
  catalog <- use_case_catalog()
  if (missing(ids) || is.null(ids)) {
    return(catalog)
  }
  dplyr::semi_join(catalog, tibble::tibble(use_case_id = ids), by = "use_case_id")
}
