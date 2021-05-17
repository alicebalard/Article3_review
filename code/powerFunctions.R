

FitBasicNoAlphaNegbin <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic without alpha")
  data$response <- data[[response]] # little trick
  HI <- data[[hybridIndex]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 A1 = paramBounds[["A1start"]],
                 Z = paramBounds[["Zstart"]])
  fit <- bbmle::mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L1, 0, HI),
                       size = SizeNegBin(A1, A1, Z, HI)),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              A1 = paramBounds[["A1LB"]],
              Z = paramBounds[["ZLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              A1 = paramBounds[["A1UB"]],
              Z = paramBounds[["ZUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

FitBasicAlphaNegbin <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model basic with alpha")
  data$response <- data[[response]] # little trick
  HI <- data[[hybridIndex]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 alpha = paramBounds[["alphaStart"]],
                 A1 = paramBounds[["A1start"]],
                 Z = paramBounds[["Zstart"]])
  fit <- bbmle::mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L1, alpha, HI),
                       size = SizeNegBin(A1, A1, Z, HI)),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              A1 = paramBounds[["A1LB"]],
              alpha = paramBounds[["alphaLB"]],
              Z = paramBounds[["ZLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              A1 = paramBounds[["A1UB"]],
              alpha = paramBounds[["alphaUB"]],
              Z = paramBounds[["ZUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

printConvergence <- function(fit) {
  convergence <- fit@details[["convergence"]]
  print(ifelse(convergence == 0, "Did converge", "Did not converge"))
}

Gtest <- function(model0, model1){
  LL0 <- Reduce(x = c(model0), f = function(accum, model){
    accum + bbmle::logLik(model)}, init = 0)
  LL1 <- Reduce(x = c(model1), f = function(accum, model){
    accum + bbmle::logLik(model)}, init = 0)
  dLL <- abs(LL1 - LL0)
  N0 <- Reduce(x = c(model0), f = function(accum, model){
    accum + length(bbmle::coef(model))}, init = 0)
  N1 <- Reduce(x = c(model1), f = function(accum, model){
    accum + length(bbmle::coef(model))}, init = 0)
  dDF <- abs(N1 - N0)
  pvalue <- 1 - stats::pchisq(2*dLL, df=dDF)
  chisqvalue <- stats::qchisq(p = pvalue, df=dDF)
  out <- data.frame(dLL = round(dLL, 2),
                    dDF = dDF,
                    pvalue = pvalue,
                    chisqvalue = chisqvalue)
  return(out)
}

getParamBounds <- function(model, data, response){
  if (model == "negbin"){
    paramBounds <- c(L1start = mean(stats::na.omit(data[[response]])),
                     L1LB = 0,
                     L1UB = max(stats::na.omit(data[[response]])),
                     L2start = mean(stats::na.omit(data[[response]])),
                     L2LB = 0,
                     L2UB = max(stats::na.omit(data[[response]])),
                     alphaStart = 0, alphaLB = -5, alphaUB = 5,
                     A1start = 10, A1LB = 1e-9, A1UB = 1000,
                     A2start = 10, A2LB = 1e-9, A2UB = 1000,
                     Zstart = 0, ZLB = -20, ZUB = 20)
  } 
  return(paramBounds)
}




myPowerAnalysis <- function(myalpha, mod, parental){
  ## Simulate data ##
  hybind = runif(n = N, min = 0, max = 1)
  simdata = data.frame(ID = 1:336, HI =hybind,
                       Sex = factor(rep(c("female", "male"), 336/2)),
                       OPG= round(parasiteLoad::MeanLoad(
                         L1 = parental, L2 = parental, alpha = myalpha,
                         hybridIndex = hybind)))
  ## Run analyis ##
  myparamBounds <- getParamBounds(model= mod, data=simdata, 
                                  response="OPG")
  
  if (mod == "negbin"){
    modAlpha <- FitBasicAlphaNegbin(data = simdata, response = "OPG", hybridIndex = "HI", 
                                    paramBounds = myparamBounds, 
                                    config =  list(optimizer = "optimx",
                                                   method = c("L-BFGS-B","bobyqa"), 
                                                   control = list(follow.on = TRUE)))
    
    modnoAlpha <- FitBasicNoAlphaNegbin(data = simdata, response = "OPG", hybridIndex = "HI", 
                                        paramBounds = myparamBounds, 
                                        config =  list(optimizer = "optimx",
                                                       method = c("L-BFGS-B","bobyqa"), 
                                                       control = list(follow.on = TRUE)))
  }
  mylrt <- Gtest(modAlpha, modnoAlpha)
  
  return(list(pval = mylrt$pvalue, hybeff = modAlpha@coef[["alpha"]]))
}

PA <- function(hybeff, whichmod, parental){
  powerAn <- data.frame(t(replicate(100, myPowerAnalysis(hybeff, whichmod, parental))))
  return(list(hybridEffectIn = hybeff, 
              power = mean(powerAn$pval < 0.05),
              hybridEffectOut = mean(unlist(powerAn$hybeff))))
}

trials <-  seq(-0.25, 0.25, 0.01)
boot_fx <- function(trial, whichmod, parental) {
  r = PA(hybeff = trial, whichmod = whichmod, parental = parental)
  res = rbind(data.frame(), r)
}

# system.time({
#   results = mclapply(trials, FUN = boot_fx, whichmod = "negbin", parental = 10, mc.cores = 20)
# })

# dfResults = data.frame(matrix(unlist(results), nrow=length(results), byrow=T))
# write.csv(dfResults, "/home/alice/Dokumente/resultsPowerAnalysisGrantEmanuel_10_negbin.csv")
