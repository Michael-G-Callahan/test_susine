if (dir.exists("../susine")) {
  options(test_susine.local_susine_path = normalizePath("../susine", winslash = "/", mustWork = TRUE))

  if ("test_susine" %in% loadedNamespaces()) {
    unloadNamespace("test_susine")
  }
  if ("susine" %in% loadedNamespaces()) {
    unloadNamespace("susine")
  }

  pkgload::load_all("../susine", quiet = TRUE, reset = TRUE)
  pkgload::load_all(".", quiet = TRUE, reset = TRUE)
}
