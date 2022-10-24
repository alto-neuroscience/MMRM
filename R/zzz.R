##################################################################
##               Utilities to run on package load               ##
##################################################################

.onLoad <- function(libname, pkgname) {
  if (requireNamespace("emmeans", quietly = TRUE)) {
    emmeans::.emm_register("mmrm", pkgname)
  }
  register_s3_method("pbkrtest", "vcovAdj", "mmrm")
  register_s3_method("nlme", "corMatrix", "corToep")
  register_s3_method("nlme", "corFactor", "corToep")
  register_s3_method("nlme", "Initialize", "corToep")
  register_s3_method("base", "coef", "corToep")


}

register_s3_method <- function(pkg, generic, class, envir = parent.frame()) {
  fun <- get(paste0(generic, ".", class), envir = envir)
  if (isNamespaceLoaded(pkg)) {
    registerS3method(generic, class, fun, envir = asNamespace(pkg))
  }
  setHook(packageEvent(pkg, "onLoad"), function(...) {
    registerS3method(generic, class, fun, envir = asNamespace(pkg))
  })
}
