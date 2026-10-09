d <- tempfile(); dir.create(d)
tar <- file.path(d, "pkg.tar.gz")
system2("R", c("CMD","build","--no-build-vignettes","ADMIT"), stdout=file.path(d,"b.log"), stderr=file.path(d,"b.log"))
f <- list.files(".", "^ADMIT_.*tar[.]gz$", full.names=TRUE)
cat("built:", basename(f), "\n")
out <- system2("R", c("CMD","check","--no-manual","--no-build-vignettes", f),
               stdout=TRUE, stderr=TRUE)
writeLines(grep("^\\*|NOTE|WARNING|ERROR", out, value=TRUE))
