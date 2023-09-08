.onLoad <- function(libname, pkgname) {

  desc <- read.dcf(file.path(system.file(package = pkgname, lib.loc = libname), "DESCRIPTION"))
  cat("Package:", desc[, "Package"], "\n")
  cat("Description:", desc[, "Description"], "\n")
  cat("Version:", desc[, "Version"], "\n")
  cat("Release Date:", desc[, "Date"], "\n")
  cat("Authors:", desc[, "Author"], "\n")
  cat("Maintainer:", desc[, "Maintainer"], "\n")
}
