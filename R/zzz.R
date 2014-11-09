.onAttach <- function(libname, pkgname) {
       version <- packageDescription("differential.coverage", field="Version")
       packageStartupMessage(paste("Welcome to differential.coverage version ", version))
}
