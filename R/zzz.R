.onAttach <- function(libname, pkgname) {
  packageStartupMessage(rep("-", 67), "\nph2rand: Design and analysis of ",
                        "randomized phase II oncology trials\n", rep("-", 67),
                        "\n\nv.0.7: For an overview of the package's ",
                        "functionality enter: ?ph2rand\n\nFor news on the ",
                        "latest updates enter: news(package = \"ph2rand\")")
}

.onUnload <- function (libpath) {
  library.dynam.unload("ph2rand", libpath)
}