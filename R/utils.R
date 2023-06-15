.onAttach <- function(libname, pkgname) {
  message("Loading package: ", pkgname)
  message("Located in library: ", libname)
  packageStartupMessage("Welcome to my package")
}
