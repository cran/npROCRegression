controlINPROCreg <-
function(
   step.p = 0.02,
   kbin = 30,
   p = 1,
   h  = c(-1,-1,-1,-1),
   seed = NULL,
   nboot = 500,
   level = 0.95,
   resample.m =c("coutcome","ncoutcome"))
   list(step.p = step.p, kbin = kbin, p = p, h = h, seed = ifelse(is.null(seed), runif(1)*1e9, seed), nboot = nboot, level=level, resample.m = match.arg(resample.m))
