.onUnload <- function (libpath) {
   library.dynam.unload("flashpcaR", libpath)
}

# https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html
.datatable.aware = TRUE

