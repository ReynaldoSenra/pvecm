#' Long-run causality in panel VECM!
#' This function applies three long-run causality tests for panel VECM proposed by Canning and Pedroni.
#'
#' @param x An object which class is 'pvecm'.
#' @param type It can be "lambda" (for the Lambda-Pearson type test), "proportion" (for the method of the proportion of units rejecting the null of no significance of the ECT) or "grouped t" (for the mean group of the ECT t-statistics method). These tests were proposed by Canning and Pedroni in a serie of papers. A description of the tests can be found in Pedroni 2019 (10.1016/b978-0-12-814367-4.00010-1) and Canning and Pedroni 2008 (10.1111/j.1467-9957.2008.01073.x).
#'
#' @import stats
#' @export
#'
#' @examples
#' data("Grunfeld", package = "plm")
#' sales <- pvecm(Grunfeld[c(3,4)], Grunfeld[5], cross_sections = Grunfeld[1], time = Grunfeld[2],
#' deterministic_long = "drift", vecm_lags = 2)
#' cause_pvecm(sales, type = "lambda")
cause_pvecm <- function(x, type = "lambda") {
  cat("
       ###                                          ###
       ###  Canning-Pedroni long-run causality test ###
       ###                                          ###\n")
  if (type == "lambda") {
    cat("\nLambda-Pearson test\n")
    chi.tests <- cbind(-2*colSums(log(x$p.value.ect)),round(pchisq(-2*colSums(log(x$p.value.ect)),2*nrow(x$p.value.ect),lower.tail = FALSE),5))
    colnames(chi.tests) <- c("Chi-statistics", "p-value")
    rownames(chi.tests) <- x$endogenous.variables
    print(chi.tests)
  } else if(type == "proportion") {
    cat("\nProportion test\n")
    proportion <- ((colSums(x$p.value.ect<0.1)*100)/nrow(x$p.value.ect)-10)/(30/sqrt(nrow(x$p.value.ect)))
    names(proportion) <- x$endogenous.variables
    print(proportion)
  } else if(type == "grouped t") {
    cat("\nGroup Mean t-test\n")
    grouped.t <- colSums(x$t.value.ect)/sqrt(nrow(x$t.value.ect))
    names(grouped.t) <- x$endogenous.variables
    print(grouped.t)
  }
}
