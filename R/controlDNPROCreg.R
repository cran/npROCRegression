controlDNPROCreg <-
function(
   step.p = 0.02,
   card.P = 50,
   link = c("probit", "logit","cloglog"),                                                         
   kbin = 30,
   p = 1,
   seed = NULL,
   nboot = 500,
   level = 0.95,                                                                
   resample.m =c("coutcome","ncoutcome"))
   list(step.p = step.p, card.P = card.P, link = match.arg(link), kbin = kbin, p = p, seed = ifelse(is.null(seed), runif(1)*1e9, seed), nboot = nboot, level = level, resample.m = match.arg(resample.m))
