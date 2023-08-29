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

.onAttach <- function(lib, pkg)  {
  packageStartupMessage("This is GoeVeg ",
                        utils::packageDescription("goeveg", field="Version"), paste(' - build: '),
                        utils::packageDate('goeveg'),
                        appendLF = TRUE)
}


.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/R-dev",
    devtools.install.args = "",
    devtools.name = "Friedemann von Lampe, Jenny Schellenberg",
    devtools.desc.author = '"Friedemann von Lampe <fvonlampe@uni-goettingen.de> [aut, cre]"',
    devtools.desc.license = "GPL (>= 2)",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])

  invisible()
}

# utils::globalVariables(c("scale_tabs"))



