#'Print a gamtest Object
#'
#'This function print the semiparametric estimation of nonlinear curves and surface.
#'@param x A gamtest object.
#'@param digits Number of digits for test statistic and p-value.
#'@param ... Other generic options.
#'@export
print.gamtest <- function (x, digits = 4, ...){
  data.bind <- x$data
  colnames(data.bind)=x$mydataname

  if (x$fcn == "gam.grptest"){
    if (ncol(x$data)==3) {
      curve.surface <- "curves"
      curve.surface2 <- "curve"
    }else {
      curve.surface <- "surfaces"
      curve.surface2 <- "surface"
    }
    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep="\n")
    cat("\n")
    cat("Comparing", x$group, "semiparametric regression", curve.surface, "\n")
    cat("Penalized semiparametric regression is used for", curve.surface2, "fitting.", "\n")
    cat("Wild-bootstrap algorithm is applied to obtain the null distribution.", "\n")
    cat("\n")
    cat("Null hypothesis: No difference between the ", x$group, " ", curve.surface, sep="", ".\n")
    cat("T = ", formatC(x$statistic, digits = digits), "   ", "p-value = ", formatC(x$p.value, digits = digits), "\n")
    cat("\n")
    invisible(x)
  }else if(x$fcn == "T.L2c"){
    if (ncol(x$data)==3) {
      curve.surface <- "curves"
      curve.surface2 <- "curve"
    }else {
      curve.surface <- "surfaces"
      curve.surface2 <- "surface"
    }
    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep="\n")
    cat("\n")
    cat("Comparing", x$group, "nonparametric regression", curve.surface, "\n")
    cat("Local polynomial regression with automatic smoothing parameter selection via", toupper(x$criterion), "is used for", curve.surface2, "fitting.", "\n")
    cat("Wild-bootstrap algorithm is applied to obtain the null distribution.", "\n")
    cat("\n")
    cat("Null hypothesis: there is no difference between the ", x$group, " ", curve.surface, sep="", ".\n")
    cat("T = ", formatC(x$statistic, digits = digits), "   ", "p-value = ", formatC(x$p.value, digits = digits), "\n")
    cat("\n")
    invisible(x)
  }else if (x$fcn == "gamm4.grptest"){
    if (ncol(x$data)==4) {
      curve.surface <- "curves"
      curve.surface2 <- "curve"
    }else {
      curve.surface <- "surfaces"
      curve.surface2 <- "surface"
    }
    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep="\n")
    cat("\n")
    cat("Comparing", x$group, "semiparametric regression", curve.surface, "\n")
    cat("Penalized semiparametric regression mixed model is used for", curve.surface2, "fitting.", "\n")
    cat("LP-bootstrap algorithm is applied to obtain the null distribution.", "\n")
    cat("\n")
    cat("Null hypothesis: there is no difference between the ", x$group, " ", curve.surface, sep="", ".\n")
    cat("T = ", formatC(x$statistic, digits = digits), "   ", "p-value = ", formatC(x$p.value, digits = digits), "\n")
    cat("\n")
    invisible(x)
  }
}
