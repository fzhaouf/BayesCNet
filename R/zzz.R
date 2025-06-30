# Package loading configuration
# This file ensures proper loading order

.onLoad <- function(libname, pkgname) {
  # Set options
  op <- options()
  op.BayesCNet <- list(
    BayesCNet.verbose = TRUE,
    BayesCNet.parallel = TRUE
  )
  toset <- !(names(op.BayesCNet) %in% names(op))
  if(any(toset)) options(op.BayesCNet[toset])

  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("BayesCNet v", utils::packageVersion("BayesCNet"))
  packageStartupMessage("Type 'citation(\"BayesCNet\")' for citing this package")
}

# Suggested collation order for DESCRIPTION file:
# Collate:
#   'AllGenerics.R'
#   'AllClasses.R'
#   'Utils.R'
#   'CreateBayesCNet.R'
#   'DataPreparation.R'
#   'RunBayesCNet.R'
#   'Accessors.R'
#   'zzz.R'
