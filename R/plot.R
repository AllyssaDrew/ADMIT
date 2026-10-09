#' Plot the efficiency curve and the two optima
#'
#' @param x An `admit_threshold` object.
#' @param show_floor Shade the region below the estimator's own noise floor, where a cut cannot be
#'   applied meaningfully. Needs `eps0` to have been passed as an [adm_eps0()] object.
#' @param ... Passed to `plot`.
#' @export
plot.admit_threshold <- function(x, show_floor = TRUE, ...) {
  tab <- x$table
  op <- par(mar = c(4.2, 4.2, 2.6, 1.2)); on.exit(par(op))
  plot(tab$cut, tab$eff, type = "n", log = "x", xlab = expression("inclusion cut on " * hat(q)),
       ylab = "NCP efficiency", main = "", panel.first = grid(col = "grey92", lty = 1), ...)

  if (show_floor && x$has_floor && length(x$floor_q)) {
    fl <- mean(x$floor_q)                              # where the estimator reads for an unadmixed person
    usr <- par("usr")
    polygon(c(10^usr[1], fl, fl, 10^usr[1]), c(usr[3], usr[3], usr[4], usr[4]),
            col = adjustcolor("grey70", 0.18), border = NA)
    text(fl, usr[3] + 0.06 * diff(usr[3:4]), " below the estimator's own floor",
         adj = 0, cex = 0.72, col = "grey35")
  }
  lines(tab$cut, tab$eff, lwd = 2.2, col = "#1F6FB4")
  d <- x$detection; s <- x$estimation
  abline(v = c(d$cut, s$cut), col = "grey45", lty = 2)
  points(d$cut, d$eff, pch = 19, col = "#1F6FB4", cex = 1.1)
  mtext(sprintf("detection %.3g   estimation %.3g   (eps0 = %.3g)", d$cut, s$cut, x$eps0),
        side = 3, line = 0.4, cex = 0.85)
  if (x$has_floor) rug(x$floor_q[x$floor_q > 0], col = adjustcolor("#C0392B", 0.25), side = 1)
  invisible(x)
}
