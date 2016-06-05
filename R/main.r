# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


### Welcome Message and Options -----

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to the GoeVeg Package")
}


.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/R-dev",
    devtools.install.args = "",
    devtools.name = "Friedemann Goral, Jenny Schellenberg",
    devtools.desc.author = '"Friedemann Goral <fgoral@gdwg.de> [aut, cre]"',
    devtools.desc.license = "GPL (>= 2)",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])

  invisible()
}



