.onAttach <- function(libname, pkgname) {
       version <- packageDescription("differential.coverage", field="Version")
       packageStartupMessage(paste("Welcome to differential.coverage version ", version)," ensemble is coming version")
}
#O Inverno Está a Сhegar 