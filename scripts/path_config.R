gart_repo_root <- function(script_file) {
  override <- Sys.getenv("GART_UC_REPO_ROOT", unset = "")
  if (nzchar(override)) return(normalizePath(override, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = FALSE)
}

gart_paths <- function(script_file) {
  repo <- gart_repo_root(script_file)
  override <- Sys.getenv("GART_UC_DATA_ROOT", unset = "")
  external <- if (nzchar(override)) normalizePath(override, winslash = "/", mustWork = FALSE) else file.path(repo, "data_external")
  list(repo_root = repo, external_data_root = external, work_root = file.path(repo, "work"), data_public_root = file.path(repo, "data_public"), results_root = file.path(repo, "results_reproduced"))
}

gart_require_file <- function(path, label = "required input") {
  if (!file.exists(path)) stop(sprintf("Required external input not found: %s (%s). See data_external/README.md or set GART_UC_DATA_ROOT.", path, label), call. = FALSE)
  path
}
