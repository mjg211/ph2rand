.onAttach <- function(libname, pkgname) {
  packageStartupMessage(rep("-", 66), "\nph2rand: Design of randomized ",
                        "comparative phase II oncology trials\n", rep("-", 66),
                        "\n\nv.0.7: For an overview of the package's ",
                        "functionality enter: ?ph2rand\n\nFor news on the ",
                        "latest updates enter: news(package = \"ph2rand\")")
}

.onUnload <- function (libpath) {
  library.dynam.unload("ph2rand", libpath)
}