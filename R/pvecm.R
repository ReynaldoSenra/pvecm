#' Panel VECM!
#' This is a function named 'pvecm' which estimate panel VECM as in Pedroni 2019 (10.1016/b978-0-12-814367-4.00010-1) and Canning and Pedroni 2008 (10.1111/j.1467-9957.2008.01073.x).
#'
#' @param right_hand_side A numeric dataframe containing all the I(1) regressors.
#' @param left_hand_side A numeric dataframe containing the dependent variable.
#' @param I_0_dols A numeric dataframe containing all the I(1) regressors. Only valid if mehtod = "D" (the long run cointegrated vector is obtained by DOLS).
#' @param cross_sections A numeric dataframe containing the index for the cross sections. It should contain consecutive integers starting from 1.
#' @param time A numeric dataframe containing the index for time.
#' @param dummies A numeric dataframe containing dummies to include in the short-run equations.
#' @param method It can be "D" if the long run cointegrated vector is obtained by DOLS or it can be "FM" if obtained by FMOLS. Default is "D". See documentation for cointReg package.
#' @param deterministic_long Deterministic variables to include in the long run cointegrated vector. It can be "none", "drift" or "trend". Default is "drift". See documentation for cointReg package.
#' @param deterministic_short Deterministic variables to include in the sort run cointegrated vector. It can be "drift" or "trend". Default is "drift". See documentation for cointReg package.
#' @param vecm_lags Lags to include in the VECM systems.
#' @param maximum_lags If you want to use the AIC, the Modified AIC or BIC information criterias, here you set the maximum number of lags that you want to consider. To compare models, you change the argument vecm_lags and use the section 'aic.each.vecm' or 'bic.each.vecm' (for example, my_vecm_estimation$aic.each.vecm) of the obtained object to compare models. Unfortunately, in the current version this need to be done manually. In other words, if you consider a maximum number of lags of 4 (maximum_lags = 4), you should set vecm_lags = 1 and after vecm_lags = 2 and vecm_lags = 3 and vecm_lags = 4 to compare models
#' @param kernel See documentation for cointReg package.
#' @param aic_small Indicates that the Modified AIC and not the AIC will be used when determined the optimal lags in the VECMS
#' @param bandwidth See documentation for cointReg package.
#' @param n.lead See documentation for cointReg package.
#' @param n.lag See documentation for cointReg package.
#' @param kmax See documentation for cointReg package.
#' @param info.crit See documentation for cointReg package.
#' @param demeaning See documentation for cointReg package. Always FALSE in the current version.
#' @param check See documentation for cointReg package.
#'
#' @import plm zoo
#' @export
#'
#' @examples
#' data("Grunfeld", package = "plm")
#' sales <- pvecm(Grunfeld[c(3,4)], Grunfeld[5], cross_sections = Grunfeld[1], time = Grunfeld[2],
#' deterministic_long = "drift", vecm_lags = 2)
#' print(sales)
pvecm <- function(right_hand_side, left_hand_side, I_0_dols = NULL, cross_sections, time, dummies = NULL, method = "D", deterministic_long = "none",
                  deterministic_short = "drift", vecm_lags = 2, maximum_lags = 3,  kernel = "ba",  aic_small = TRUE,
                  bandwidth = "and", n.lead = NULL, n.lag = NULL,
                  kmax = "k4", info.crit = "AIC", demeaning = FALSE,
                  check = TRUE) {

  vcem_cross <- list()
  individual_cross <- list()
  resid_cd <- list()
  pmg_model <- list()
  vecm_output <- list()
  ECT <- c()
  data <- base::cbind(cross_sections,time,left_hand_side, right_hand_side)
  if (!is.null(dummies)) {
    data <- base::cbind(cross_sections,time,left_hand_side,right_hand_side,dummies)
  }
  if (!is.null(I_0_dols)) {
    data <- base::cbind(cross_sections,time,left_hand_side, I_0_dols, right_hand_side)
    if (!is.null(dummies)) {
      data <- base::cbind(cross_sections,time,left_hand_side, I_0_dols, right_hand_side,dummies)
    }
  }
  if (any(is.na(data) | is.infinite(base::as.matrix(data)))) stop('There are NAs or Infinity in the data')
  dato <- base::as.data.frame(data)
  cr <- dplyr::count(dplyr::distinct(dato, dato[1]))[[1]]
  for (i in 1:cr){
    right_hand_sides <- base::subset(dato[c(base::colnames(cross_sections), base::colnames(right_hand_side))], dato[[1]] %in% i)
    right_hand_sides <- right_hand_sides[-1]
    left_hand_sides <-  base::subset(dato[c(base::colnames(cross_sections), base::colnames(left_hand_side))], dato[[1]] %in% i)
    left_hand_sides <-  left_hand_sides[-1]
    datas <- base::diff(base::as.matrix(base::cbind(left_hand_sides, right_hand_sides)))
    if (method == "D") {
      drift <- rep(1, times = nrow(right_hand_sides))
      I_0_dolss <-  I_0_dols
      if (!is.null(I_0_dols)) {
        I_0_dolss <-  base::subset(dato[c(base::colnames(cross_sections), base::colnames(I_0_dols))], dato[[1]] %in% i)
        I_0_dolss <-  I_0_dolss[-1]
        datas <- base::diff(base::as.matrix(base::cbind(left_hand_sides, I_0_dolss, right_hand_sides)))
        drift <- base::cbind(drift, I_0_dolss)
      }
      if (deterministic_long == "drift") {
        coefficient_names <- c("drift",base::colnames(datas)[c(2:base::ncol(datas))])
        long_run <- cointReg::cointReg(right_hand_sides, left_hand_sides, deter = drift, method = "D", kernel = kernel,
                                       bandwidth = bandwidth, n.lead = n.lead, n.lag = n.lag, kmax = kmax, info.crit = info.crit, demeaning = demeaning, check = check)
        if (i == 1) {
          matrix_coef <- matrix(,nrow = cr,ncol = length(long_run$theta))
          matrix_t <- matrix(,nrow = cr,ncol = length(long_run$t.theta))
        }
        matrix_coef[i,] <- long_run$theta
        matrix_t[i,] <- long_run$t.theta
        ECT <- append(ECT,long_run$residuals)
      } else if (deterministic_long == "trend") {
        coefficient_names <- c("trend", "drift", base::colnames(datas)[c(2:base::ncol(datas))])
        trend <- c(1:length(drift))
        deter <- base::cbind(trend,drift)
        long_run <- cointReg::cointReg(right_hand_sides, left_hand_sides, deter = deter, method = "D", kernel = kernel,
                                       bandwidth = bandwidth, n.lead = n.lead, n.lag = n.lag, kmax = kmax, info.crit = info.crit, demeaning = demeaning, check = check)
        if (i == 1) {
          matrix_coef <- matrix(,nrow = cr,ncol = length(long_run$theta))
          matrix_t <- matrix(,nrow = cr,ncol = length(long_run$t.theta))
        }
        matrix_coef[i,] <- long_run$theta
        matrix_t[i,] <- long_run$t.theta
        ECT <- append(ECT,long_run$residuals)
      } else {
        coefficient_names <- base::colnames(datas)[c(2:base::ncol(datas))]
        long_run <- cointReg::cointReg(right_hand_sides, left_hand_sides, deter = I_0_dolss, method = "D", kernel = kernel,
                                       bandwidth = bandwidth, n.lead = n.lead, n.lag = n.lag, kmax = kmax, info.crit = info.crit, demeaning = demeaning, check = check)
        if (i == 1) {
          matrix_coef <- matrix(,nrow = cr,ncol = length(long_run$theta))
          matrix_t <- matrix(,nrow = cr,ncol = length(long_run$t.theta))
        }
        matrix_coef[i,] <- long_run$theta
        matrix_t[i,] <- long_run$t.theta
        ECT <- append(ECT,long_run$residuals)
      }
    } else {
      if (deterministic_long == "drift") {
        deter <- rep(1, times = nrow(right_hand_sides))
        long_run <- cointReg::cointReg(right_hand_sides, left_hand_sides, deter = deter, method = "FM", kernel = kernel,
                                       bandwidth = bandwidth, demeaning = demeaning, check = check)
        coefficient_names <- c("drift",base::colnames(datas)[c(2:base::ncol(datas))])
        if (i == 1) {
          matrix_coef <- matrix(,nrow = cr,ncol = length(long_run$theta))
          matrix_t <- matrix(,nrow = cr,ncol = length(long_run$t.theta))
        }
        matrix_coef[i,] <- long_run$theta
        matrix_t[i,] <- long_run$t.theta
        ECT <- append(ECT,long_run$residuals)
      } else if (deterministic_long == "trend") {
        coefficient_names <- c("drift","trend", base::colnames(datas)[c(2:base::ncol(datas))])
        drift <- rep(1, times = nrow(right_hand_sides))
        trend <- c(1:length(drift))
        deter <- base::cbind(drift, trend)
        long_run <- cointReg::cointReg(right_hand_sides, left_hand_sides, deter = deter, method = "FM", kernel = kernel,
                                       bandwidth = bandwidth, demeaning = demeaning, check = check)
        if (i == 1) {
          matrix_coef <- matrix(,nrow = cr,ncol = length(long_run$theta))
          matrix_t <- matrix(,nrow = cr,ncol = length(long_run$t.theta))
        }
        matrix_coef[i,] <- long_run$theta
        matrix_t[i,] <- long_run$t.theta
        ECT <- append(ECT,long_run$residuals)
      } else {
        coefficient_names <- base::colnames(datas)[c(2:base::ncol(datas))]
        long_run <- cointReg::cointReg(right_hand_sides, left_hand_sides, method = "FM", kernel = kernel,
                                       bandwidth = bandwidth, demeaning = demeaning, check = check)
        if (i == 1) {
          matrix_coef <- matrix(,nrow = cr,ncol = length(long_run$theta))
          matrix_t <- matrix(,nrow = cr,ncol = length(long_run$t.theta))
        }
        matrix_coef[i,] <- long_run$theta
        matrix_t[i,] <- long_run$t.theta
        ECT <- append(ECT,long_run$residuals)
      }
    }
    # ECT1 <- lag(long_run$residuals)
    # ECT1 <- ECT1[c(2:length(ECT1))]
    ECT1 <- as.data.frame(utils::head(long_run$residuals,-1))  # changed for the package
    colnames(ECT1) <- "ECT1"
    if (i == 1) {
      ecm_t <- matrix(,nrow = cr,ncol = base::ncol(datas))
      ecm_pvalue <- matrix(,nrow = cr,ncol = base::ncol(datas))
      aic <- matrix(,nrow = cr,ncol = base::ncol(datas))     # new for the package
      bic <- matrix(,nrow = cr,ncol = base::ncol(datas))     # new for the package
      observat <- matrix(,nrow = cr,ncol = base::ncol(datas))   # new for the package
    }
    if (!is.null(dummies)) {
      dummys <-  base::subset(dato[c(base::colnames(cross_sections), base::colnames(dummies))], dato[[1]] %in% i)
      dummys <-  dummys[-1]
      dummys <- base::as.data.frame(dummys[-1,])
      base::colnames(dummys) <- base::colnames(dummies)
      datas <- base::cbind(datas,dummys,ECT1)
      if (vecm_lags != maximum_lags) {
        datus <- stats::ts(datas[c(-(maximum_lags-vecm_lags):-1),])  # new for the package
      } else {
        datus <- stats::ts(datas)  # new for the package
      }
      datas <- stats::ts(datas)
      names_vecm <- base::colnames(datas)
      if (i == 1) {
        ecm_t <- matrix(,nrow = cr,ncol = (base::ncol(datas)-base::ncol(dummies)-1))
        ecm_pvalue <- matrix(,nrow = cr,ncol = (base::ncol(datas)-base::ncol(dummies)-1))
        aic <- matrix(,nrow = cr,ncol = (base::ncol(datas)-base::ncol(dummies)-1))    # new for the package
        bic <- matrix(,nrow = cr,ncol = (base::ncol(datas)-base::ncol(dummies)-1))    # new for the package
        observat <- matrix(,nrow = cr,ncol = (base::ncol(datas)-base::ncol(dummies)-1))    # new for the package
      }
      for (j in 1:(base::ncol(datas)-base::ncol(dummies)-1)) {
        if (deterministic_short == "drift") {
          formulation <- stats::as.formula(paste(names_vecm[j],paste(paste(paste("L(", names_vecm[-c((length(names_vecm)-base::ncol(dummies)):length(names_vecm))],",1:",vecm_lags,")", sep = ""), collapse = "+"), paste(names_vecm[c((length(names_vecm)-base::ncol(dummies)):length(names_vecm))],collapse = "+")  , sep = "+" ), sep = "~"))

          short_run <- dynlm::dynlm(formulation, data = datas)
          short_run1 <- dynlm::dynlm(formulation, data = datus)  # new for the package

        } else {

          short_run <- dynlm::dynlm(stats::as.formula(paste(names_vecm[j],paste("trend(datas[,1])",paste(paste("L(", names_vecm[-c((length(names_vecm)-base::ncol(dummies)):length(names_vecm))],",1:",vecm_lags,")", sep = ""), collapse = "+"),paste(names_vecm[(length(names_vecm)-base::ncol(dummies)):length(names_vecm)],collapse = "+"), sep = "+" ) ,sep = "~")), data = datas)
          short_run1 <- dynlm::dynlm(stats::as.formula(paste(names_vecm[j],paste("trend(datas[,1])" ,paste(paste("L(", names_vecm[-c((length(names_vecm)-base::ncol(dummies)):length(names_vecm))],",1:",vecm_lags,")", sep = ""), collapse = "+"),paste(names_vecm[(length(names_vecm)-base::ncol(dummies)):length(names_vecm)],collapse = "+"), sep = "+" ) ,sep = "~")), data = datus)  # new for the package

        }
        ecm_t[i,j]<-summary(short_run)$coefficients[nrow(summary(short_run)$coefficients),3]
        ecm_pvalue[i,j]<-summary(short_run)$coefficients[nrow(summary(short_run)$coefficients),4]
        aic[i,j] <- AICcmodavg::AICc(short_run1,second.ord = aic_small)    # new for the package
        bic[i,j] <- BIC(short_run1)    # new for the package
        observat[i,j] <- stats::nobs(short_run1)   # new for the package
        vcem_cross[[j]] <- list(names_vecm[j], summary(short_run)$coefficients)
      }
    } else {
      datas <- base::cbind(datas,ECT1)
      if (vecm_lags != maximum_lags) {
        datus <- stats::ts(datas[c(-(maximum_lags-vecm_lags):-1),])  # new for the package
      } else {
        datus <- stats::ts(datas)  # new for the package
      }
      datas <- stats::ts(datas)
      names_vecm <- base::colnames(datas)
      for (j in 1:(base::ncol(datas)-1)) {
        if (deterministic_short == "drift") {
          formulation <- stats::as.formula(paste(names_vecm[j],paste(paste(paste("L(", names_vecm[-length(names_vecm)],",1:",vecm_lags,")", sep = ""), collapse = "+"),names_vecm[length(names_vecm)], sep = "+" ) ,sep = "~"))

          short_run <- dynlm::dynlm(formulation, data = datas)
          short_run1 <- dynlm::dynlm(formulation, data = datus)  # new for the package
        } else {

          short_run <- dynlm::dynlm(stats::as.formula(paste(names_vecm[j],paste("trend(datas[,1])" ,paste(paste("L(", names_vecm[-length(names_vecm)],",1:",vecm_lags,")", sep = ""), collapse = "+"),names_vecm[length(names_vecm)], sep = "+" ) ,sep = "~")), data = datas)
          short_run1 <- dynlm::dynlm(stats::as.formula(paste(names_vecm[j],paste("trend(datas[,1])" ,paste(paste("L(", names_vecm[-length(names_vecm)],",1:",vecm_lags,")", sep = ""), collapse = "+"),names_vecm[length(names_vecm)], sep = "+" ) ,sep = "~")), data = datus)
        }
        ecm_t[i,j]<-summary(short_run)$coefficients[nrow(summary(short_run)$coefficients),3]
        ecm_pvalue[i,j]<-summary(short_run)$coefficients[nrow(summary(short_run)$coefficients),4]
        aic[i,j] <- AICcmodavg::AICc(short_run1,second.ord = aic_small)     # new for the package
        bic[i,j] <- BIC(short_run1)     # new for the package
        observat[i,j] <- stats::nobs(short_run1)    # new for the package
        vcem_cross[[j]] <- summary(short_run)$coefficients
        names(vcem_cross)[[j]] <- names_vecm[j]
      }
    }
    individual_cross[[i]]  <- vcem_cross
  }
  stata_mg <- plm::pdata.frame(base::cbind(data,ECT), index = c(base::colnames(cross_sections),base::colnames(time)))
  data_pmg <- plm::pdata.frame(base::cbind(data,ECT), index = c(base::colnames(cross_sections),base::colnames(time)), drop.index = TRUE)
  names_vecm_pmg <-names(data_pmg)
  cc <- vecm_lags
  names_pmg <- c()
  if (!is.null(dummies)) {
    for (i in 1:vecm_lags) {
      names_pmg[i] <- paste(paste("lag(diff(", names_vecm_pmg[-c((length(names_vecm_pmg)-base::ncol(dummies)):length(names_vecm_pmg))],  paste("),",toString(i),"L)", sep = ""), sep = ""), collapse = "+")
    }
    right_hand <- paste(names_pmg, collapse = "+")
    dummes <- paste(names(dummies),collapse = "+")
    for (j in 1:(base::ncol(datas)-base::ncol(dummies)-1)) {
      formulations <- stats::as.formula(paste(paste("diff(",names_vecm_pmg[j],")", sep = ""), paste(right_hand, dummes, paste("lag(",names_vecm_pmg[length(names_vecm_pmg)],")", sep = ""),sep = "+"), sep = "~"))
      if (deterministic_short == "trend") {
        equation_pmg <- pmg(formulations, trend = TRUE, data = data_pmg)
        pmg_model[[j]] <- summary(equation_pmg)
        resid_cd[[j]] <- pcdtest(equation_pmg, test = "cd")
      } else {
        equation_pmg <- pmg(formulations, data = data_pmg)
        pmg_model[[j]] <- summary(equation_pmg)
        resid_cd[[j]] <- pcdtest(equation_pmg, test = "cd")
      }
    }
  } else {
    for (i in 1:vecm_lags) {
      names_pmg[i] <- paste(paste("lag(diff(", names_vecm_pmg[-length(names_vecm_pmg)],  paste("),",toString(i),"L)", sep = ""), sep = ""), collapse = "+")
    }
    right_hand <- paste(names_pmg, collapse = "+")
    for (j in 1:(base::ncol(datas)-1)) {
      formulations <- stats::as.formula(paste(paste("diff(",names_vecm_pmg[j],")", sep = ""), paste(right_hand,  paste("lag(",names_vecm_pmg[length(names_vecm_pmg)],")", sep = ""),sep = "+"), sep = "~"))
      if (deterministic_short == "trend") {
        equation_pmg <- pmg(formulations, trend = TRUE, data = data_pmg)
        pmg_model[[j]] <- summary(equation_pmg)
        resid_cd[[j]] <- pcdtest(equation_pmg, test = "cd")
      } else {
        equation_pmg <- pmg(formulations, data = data_pmg)
        pmg_model[[j]] <- summary(equation_pmg)
        resid_cd[[j]] <- pcdtest(equation_pmg, test = "cd")
      }
    }
  }
  names(pmg_model) <- names_vecm_pmg[c(1:length(pmg_model))]
  long_run_coeff <- base::cbind(colMeans(matrix_coef),colSums(matrix_t)/sqrt(nrow(matrix_t)))
  base::colnames(long_run_coeff) <- c("coefficient", "t-statistics")
  rownames(long_run_coeff) <- coefficient_names
  vecm_output <- list(estimation.method = method, cross.section = cr, endogenous.variables = names_vecm_pmg[c(1:length(pmg_model))], lags.vecm = vecm_lags, individual.coefficients = matrix_coef, individual.coefficients.t.statistics = matrix_t, individual.vecm = individual_cross, long.run.vector = long_run_coeff, VECMS = pmg_model, observations.each.vecm = observat, aic.each.vecm = aic, bic.each.vecm = bic, p.value.ect = ecm_pvalue, t.value.ect = ecm_t, resid.cd.test = resid_cd, first.step.residuals = ECT)
    class(vecm_output) <- c("pvecm", class(vecm_output))
  return(vecm_output)
}
