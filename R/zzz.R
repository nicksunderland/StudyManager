.onLoad <- function(libname, pkgname){

  gwama_exe <- system.file("bin", "GWAMA_v2", "GWAMA", package=pkgname)

  if(gwama_exe == "") {

    rlog::log_warn("No GWAMA execute able found, attempting to make...")

    gwama_dir <- system.file("bin", "GWAMA_v2", package=pkgname)

    cmd <- paste0("make --directory=", gwama_dir)

    system(cmd)

  }

  gwama_exe <- system.file("bin", "GWAMA_v2", "GWAMA", package=pkgname)

  if(gwama_exe != "") {

    rlog::log_info(paste("GWAMA executable located:", gwama_exe))

  }

}
