.onAttach <- function(libname, pkgname) {
  packageStartupMessage("  ", rep("-", 83), "\n  ph2rand: Design of randomized",
                        " comparative phase II oncology trials with a ",
                        "Bernoulli\n           distributed primary outcome\n  ",
                        rep("-", 83), "\n\n  For an overview of the package's ",
                        "functionality enter: ?ph2rand\n\n  For news on the ",
                        "latest updates enter: news(package = \"ph2rand\")")
}

.onUnload <- function (libpath) {
  library.dynam.unload("ph2rand", libpath)
}