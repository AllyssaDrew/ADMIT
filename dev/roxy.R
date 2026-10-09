if (!requireNamespace("roxygen2", quietly = TRUE)) { cat("roxygen2: NOT INSTALLED\n"); quit(status=0) }
cat("roxygen2", as.character(packageVersion("roxygen2")), "\n")
roxygen2::roxygenise("ADMIT", roclets = c("rd", "namespace"))
cat("man pages now:", length(list.files("ADMIT/man")), "\n")
cat(paste(" ", list.files("ADMIT/man")), sep="\n")
