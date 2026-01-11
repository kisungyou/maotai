## RETICULATE : global reference
.onLoad <- function(libname, pkgname) {
  if (!interactive()) {
    # Tell rgl to use an off-screen "null" device in headless checks
    options(rgl.useNULL = TRUE)
    Sys.setenv(RGL_USE_NULL = "TRUE")
  }
}

# .pkgenv <- new.env(parent = emptyenv())
.onAttach <- function(...) {
  if (interactive()) {
    packageStartupMessage("* maotai ", packageVersion("maotai"),
                          " — Tools for Matrix Algebra, Optimization and Inference")
  }
}
# .onAttach <- function(...){
#   ## Retrieve Year Information
#   date <- date()
#   x <- regexpr("[0-9]{4}", date)
#   this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
#   
#   # Retrieve Current Version
#   this.version = packageVersion("maotai")
#   
#   ## Print on Screen
#   packageStartupMessage("**-----------------------------------------------------------------**")
#   packageStartupMessage("** maotai")
#   packageStartupMessage("**  - Tools for Matrix Algebra, Optimization and Inference Problems")
#   packageStartupMessage("**")
#   packageStartupMessage("** Version    : ",this.version,"      (",this.year,")",sep="")
#   packageStartupMessage("** Maintainer : Kisung You (kisung.you@outlook.com)")
#   packageStartupMessage("** Website    : https://www.kisungyou.com/maotai")
#   packageStartupMessage("**")
#   packageStartupMessage("** Please share any bugs or suggestions to the maintainer.")
#   packageStartupMessage("**-----------------------------------------------------------------**")
# }

.onUnload <- function(libpath) {
  library.dynam.unload("maotai", libpath)
}
