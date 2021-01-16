
Break = "\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = "   Welcome to 'EDP 380C.16 Hierarchical Linear Modeling'.
   Programs developed by Reza Norouzian, Copyright (C) 2019-present"

message(Break, notice, Break)


#============================================================================================================================

#palette('R3')

#=============================================================================================================================

rm.allrowNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[rowSums(is.na(i) | i == "") != ncol(i), , drop = FALSE])
    
  } else { X[rowSums(is.na(X) | X == "") != ncol(X), , drop = FALSE] }
}

#=============================================================================================================================           

rm.allcolNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[, colSums(is.na(i) | i == "") != nrow(i), drop = FALSE])
    
  } else { X[, colSums(is.na(X) | X == "") != nrow(X), drop = FALSE] }
}           
      
#=============================================================================================================================
           
rm.colrowNA <- function(X){

r <- rm.allrowNA(X)
rm.allcolNA(r)  

}           
                  
#=============================================================================================================================           
           

trim <- function(X){
  X <- setNames(X, trimws(names(X)))
  y <- sapply(names(X), function(x) is.character(as.vector(X[[x]])))
  X[y] <- lapply(X[y], trimws)
  return(X)
}

#=============================================================================================================================

data.str <- function(dat, drop = NULL){
  
  dat <- if(is.null(drop)) dat else drop.col(dat, vec = drop)
  
  setNames(lapply(names(dat), function(i) sort(unique(dat[[i]]))), names(dat))
}                                          

#=============================================================================================================================

long <- function(data, one.cols, multi.cols = NULL, multi.colnames = NULL, time.name = "time", 
                 time.val = NULL, replace = "#N/A", with = NA, drop = NULL, order.by = one.cols[1]){
  
  data <- rm.colrowNA(trim(data)) 
  
  data[sapply(data, `%in%`, replace)] <- with
  
  varying <- if(is.null(multi.cols)) setdiff(names(data), one.cols) else multi.cols
  
  time.val <- if(is.null(time.val) & !is.null(multi.cols)) seq_along(varying[[1L]]) else time.val
  
  if(is.null(time.val) || length(time.val) == 1) stop("Please provide unique values for 'time.val'.", call. = FALSE)
  
  res <- tryCatch({na.omit(reshape(data, direction = 'long', timevar = time.name, times = time.val,
                                   idvar = one.cols, v.names = multi.colnames,
                                   varying = varying, drop = drop))}, 
                  error = function(e) {
                    if(as.character(e) == "Error in guess(varying): failed to guess time-varying variables from their names\n") stop("Please provide the 'multi.cols'.", call. = FALSE) 
                    else e
                  })
  
  if(inherits(res, "simpleError")) stop(paste0(res[[2]]), call. = FALSE)
  res <- res[order(res[[order.by]]), ]
  
  rownames(res) <- NULL
  
  res[] <- lapply(res, function(x) type.convert(as.character(x), as.is = TRUE))
  return(res)
}                

#=============================================================================================================================

mask <- function(data, what, full = FALSE){
  
  data[] <- lapply(data, function(x) type.convert(as.character(x), as.is = TRUE))
  
  f1 <- function(x) as.numeric(factor(x, levels = unique(x)))
  f2 <- function(x) {
    temp <- substr(x, 1, 1)
    paste0(temp, ave(x, temp, FUN = function(y) match(y, unique(y))))
  }
  cols <- names(data)[sapply(data, is.numeric)]
  num.cols <- cols[cols %in% what]
  cols <- names(data)[sapply(data, is.character)]
  char.cols <- cols[cols %in% what]
  
  if(!full){
    
    if(length(num.cols))  data[num.cols] <- lapply(data[num.cols], f1)
    if(length(char.cols)) data[char.cols] <- lapply(data[char.cols], f2)
    
  }else{
    
    data[what] <- lapply(data[what], f1)
  }
  return(data)
}

#=============================================================================================================================

make.dummy <- function(data, what){
  
  ff <- function(data, what){
    
    if(!inherits(data, 'data.frame')) stop("'data' must be a data.frame.", call. = FALSE)
    if(!is.character(what)) stop("'what' must be a character.", call. = FALSE)
    what <- trimws(what)
    data <- trim(data)
    if(!is.character(as.vector(data[[what]]))) stop('Not a character variable.', call. = FALSE)
    formula <- as.formula(paste('~', what))
    data.frame(model.matrix(formula, data = data))[-1L]
  }
  
  cbind(data, do.call(cbind, lapply(what, ff, data = data)))
}                     

#=============================================================================================================================
                     
rm.terms <- function(form, term) {
    fterms <- terms(form)
    fac <- attr(fterms, "factors")
    idx <- which(as.logical(fac[term, ]))
    new.fterms <- stats::drop.terms(fterms, dropx = idx, keep.response = TRUE)
    return(as.formula(new.fterms))
  }                     
                     
#=============================================================================================================================

post.mixed <- function(fit, formula = NULL, plot = TRUE, by = NULL, var = NULL, horiz = TRUE, adjust = "tukey", type = "response", compare = FALSE, weights = c("equal", "proportional", 
                                                                                                                                                                "outer", "cells", "flat", 
                                                                                                                                                                "show.levels")[1], ...){
  
  limit <- nobs(fit)
  tm <- terms(fit)
  
  f <- if(is.null(formula)) as.formula(bquote(pairwise ~ .(tm[[3L]]))) else as.formula(formula)
  
  av <- emmeans::.all.vars(tm)[-1L]
  
  cl <- if(inherits(fit, "lme")) attr(tm, "dataClasses")[-1L] else setNames(sapply(av, function(i) class(model.frame(fit)[[i]])), av)
  
  all.factor <- all(cl == "factor") || all(cl == "character") || all(cl == "character" | cl == "factor")
  
  no.contrast <- length(av) == 1L & !all.factor
  
  ems <- if(all.factor || no.contrast)  { eval(substitute(emmeans::emmeans(fit, f, infer = c(TRUE, TRUE), type = type, pbkrtest.limit = limit, weights = weights))) 
    
  } else {
    
    var <- if(is.null(var)) names(sort(cl[grep("integer|numeric", cl)])[1]) else var
    
    f <- rm.terms(f, var)
    
    eval(substitute(emmeans::emtrends(fit, f, var = var, infer = c(TRUE, TRUE), type = type, pbkrtest.limit = limit, weights = weights)))
  } 
  
  xlab <- if(!all.factor)  paste(var, "(slope)") else "Estimated Means"
  
  print(plot(ems, by = by, comparisons = compare, horizontal = horiz, adjust = adjust, xlab = xlab, ...))
  
  em <- as.data.frame(ems[[2]])
  
  out <- if(inherits(fit, c("lmerMod", "lmerModLmerTest", "lme4", "lme")) & !no.contrast){
    
    vc <- VarCorr(fit)
    
    sigma <- if(inherits(fit, "lme")) sqrt(sum(as.numeric(vc[,"Variance"]), na.rm = TRUE)) else sqrt(sum(as.numeric(c(attr(vc[[1]], "stddev"), attr(vc, "sc")))^2, na.rm = TRUE))
    
    edf <- min(as.data.frame(ems[[1L]])$df, na.rm = TRUE)
    
    ef <- as.data.frame(emmeans::eff_size(ems[[1L]], sigma = sigma, edf = edf))[c(2,5,6)]
    
    res <- cbind(em, ef)
    names(res)[c(2,5:7, 9:11)] <- c(if(all.factor) "mean.dif" else paste0("slope.dif.", var), "lower", "upper", "t.value", "Cohen.d", "lower.d", "upper.d")
    list(emmeans = ems[[1]], contrasts = roundi(res, 3))
    
  }
  
  else {
    
    ems
  }
  
  return(out)
}      

#=================================================================================================================================
                                                                       
lo_ave_up <- function(data = NULL, vars, vals = NULL){
  
  if(is.null(vals)){
    sapply(vars, function(x) 
      round(setNames(mean(data[[x]]) + c(-1, 0, 1)*sd(data[[x]]), 
               paste0(x, c('-1SD', '.Mean', '+1SD'))), 9), simplify = FALSE) 
  } else {
    
    setNames(lapply(vars, function(i) vals), vars)
  }
}                                    
                                                                        
#=================================================================================================================================
                                  
data.tree <- function(data, toplab = NULL, cex = 1, auto = FALSE, ...){
  
  if(auto){    
    cats <- sapply(data, Negate(is.numeric))  
    data <- data[cats]
  }
  
  toplab <- if(is.null(toplab)) names(data) else toplab
  
  plotrix::sizetree(data, toplab = toplab, stacklabels = FALSE, border = 0, base.cex = cex, ...)
  
}                                 
                                  
#=================================================================================================================================
  
plot.prof <- function(fit){
  
  if(!inherits(fit, c("lmerMod", "lmerModLmerTest", "lme4", "glmmTMB", "glmerMod"))) stop("Model not supported.", call. = FALSE)
  
  pp <- profile(fit, signames = FALSE)
  
  dd <- as.data.frame(pp)
  
  if(".zeta" %in% names(dd)) names(dd)[which(names(dd) == ".zeta")] <- "value"
  if(inherits(fit, "glmmTMB")) dd$value <- sqrt(dd$value)
  
  ggplot2::ggplot(dd,aes(.focal, value)) +  geom_hline(yintercept = 0, colour = 8, linetype = 2) +
    geom_line(colour = 2) + geom_point(colour = 2) +
    facet_wrap(~.par, scale = "free_x") + xlab("Parameter Value") +
    ylab("Zeta (Normal)")
}

                                    
#=============================================================================================================================
      
                                    
latent.reg <- function(fit, formula, group.id, std = FALSE,
                       sing.check = TRUE, ridge.constant = 0, tol = 1e-12) 
{
  R <- SDs <- matrix()
  if (regexpr("merMod", class(fit)) > 0) {
    R <- attr(lme4::VarCorr(fit)[[group.id]], "correlation")
    SDs <- diag(attr(lme4::VarCorr(fit)[[group.id]], "stddev"))
  } else if (class(fit) == "lme") {
    R <- nlme::corMatrix(fit$modelStruct[[1]])[[group.id]]
    SDs <- diag(fit$sigma * attr(nlme::corMatrix(fit$modelStruct[[1]])[[group.id]], "stdDev"))
  } else stop("fitted model must be lme4 or nlme fit", call. = FALSE)
  vars <- colnames(R)
  formula <- as.character(formula)
  out_loc <- which(vars == formula[2])
  pred.names <- trimws(strsplit(formula[3], "+", fixed = TRUE)[[1]])
  pred_loc <- which(vars %in% pred.names)
  R <- stats::cov2cor(R + ridge.constant * diag(nrow(R)))
  if (tol < 1e-12) warning("Do not reduce tol, solution may not be numerically stable.")
  if (sing.check) {
    if (Matrix::det(R) < tol) {
      stop("Singular covariance matrix, did 'lme4' converge?", call. = FALSE)
    }
  } else {
    warning("Solution may not be numerically stable.", call. = FALSE)
  }
  res <- array(dim = length(pred_loc))
  if (std) {
    res <- base::solve(R[pred_loc, pred_loc], R[out_loc, pred_loc], tol = tol)
  } else {
    C <- SDs %*% R %*% SDs
    res <- base::solve(C[pred_loc, pred_loc], C[out_loc, pred_loc], tol = tol)
  }
  names(res) <- if(length(pred.names) > 1) colnames(R[pred_loc, pred_loc]) else pred.names
  return(res)
}

#=================================================================================================================================================


latent.me <- function (fit, formula, group.id, std = TRUE, digits = 3, prog.bar = "none",
                       sing.check = TRUE, ridge.constant = 1e-06, tol = 1e-12,
                       seed = 123, nsim = 499, level = .95, parallel = c("no", "multicore", "snow")[2], ncpus = 7) {
  
  formals(latent.reg)[formalArgs(latent.reg)[-1L]] <- c(
    formula, group.id, std, sing.check, ridge.constant, tol)
  
  boot.res <- list()
  if (regexpr("merMod", class(fit)) > 0) {
    boot.res <- lme4::bootMer(
      fit, nsim = nsim, seed = seed, .progress = prog.bar,
      parallel = parallel, ncpus = ncpus,
      FUN = latent.reg)
  } else if (class(fit) == "lme") {
    boot.res <- lmeresampler::parametric_bootstrap(fit, latent.reg, nsim)
  } else stop("fitted model must be lme4 or nlme fit", call. = FALSE)
  
  se.s <- apply(boot.res$t, 2, sd, na.rm = TRUE)
  z.value <- abs((2 * boot.res$t0 - colMeans(boot.res$t, na.rm = TRUE)) / se.s)
  
  boot.p <- 2 * pnorm(-z.value)
  
  lower <- upper <- array()
  if (regexpr("merMod", class(fit)) > 0) {
    ci <- confint(boot.res, level = level, type = "norm")
    lower <- ci[, 1]
    upper <- ci[, 2]
  } else if (class(fit) == "lme") {
    ci <- sapply(1:length(se.s),
                 function (i) boot::boot.ci(boot.res, conf = level, type = "norm", index = i)$normal)
    lower <- ci[2, ]
    upper <- ci[3, ]
  }
  
  round(data.frame(
    est. = boot.res$t0, SE = se.s, z.value = z.value, p.value = boot.p,
    lower = lower, upper = upper), digits)
}        

#=================================================================================================================================
                                    
cor2cov <- function (R, sds, names = NULL) 
{
    p <- (d <- dim(R))[1L]
    if (!is.numeric(R) || length(d) != 2L || p != d[2L]) 
        stop("'V' is not a square numeric matrix")
    if (any(!is.finite(sds))) 
        warning("sds had 0 or NA entries; non-finite result is doubtful")
    if (p != length(sds)) 
        stop("The standard deviation vector and correlation matrix have a different number of variables")
    S <- R
    S[] <- sds * R * rep(sds, each = p)
    if (!is.null(names)) {
        stopifnot(length(names) == p)
        rownames(S) <- colnames(S) <- names
    }
    S
}
                                    
#================================================================================================================================
                                    
n_par <- function(...){
  
  f <- function(fit){
    total <- attr(logLik(fit), "df")
    fixed <- length(fixef(fit))
    random <- ncol(ranef(fit))
    error <- if(inherits(fit, "lme") || isLMM(fit)) 1 else 0
    variance <- if(inherits(fit, "lme")) length(coef(fit$modelStruct$varStruct, unconstrained=FALSE)) else 0
    res.cov <- if(inherits(fit, "lme")) length(coef(fit$modelStruct$corStruct, unconstrained=FALSE)) else 0
    return(c(total = total, fixed = fixed, random = random, error = error, variance = variance,
             res.cov = res.cov))
  }
  
  m <- as.data.frame(purrr::map_dfr(.x = list(...), .f = f))
  rownames(m) <- substitute(...())
  m
}                                         
 
#=================================================================================================================================

# Some hack to turn off unnneeded tick mark on the 3rd and 4th axes of plot effects
                                    
plot.efflist <- function (x, selection, rows, cols, graphics = TRUE, 
                          lattice, ...) 
{
  lattice <- if (missing(lattice)) 
    list()
  else lattice
  if (!missing(selection)) {
    if (is.character(selection)) 
      selection <- gsub(" ", "", selection)
    return(plot(x[[selection]], lattice = lattice, ...))
  }
  effects <- gsub(":", "*", names(x))
  neffects <- length(x)
  mfrow <- mfrow(neffects)
  if (missing(rows) || missing(cols)) {
    rows <- mfrow[1]
    cols <- mfrow[2]
  }
  for (i in 1:rows) {
    for (j in 1:cols) {
      if ((i - 1) * cols + j > neffects) 
        break
      more <- !((i - 1) * cols + j == neffects)
      lattice[["array"]] <- list(row = i, col = j, 
                                 nrow = rows, ncol = cols, more = more)
      pp <- plot(x[[(i - 1) * cols + j]], lattice = lattice, 
                 ...)
      # hack to turn off opposite side tick marks
      pp$x.scales$tck=c(1,0)
      pp$y.scales$tck=c(1,0)
      print(pp)
    }
  }
}
environment(plot.efflist) <- asNamespace("effects")
                                    
#=================================================================================================================================
                                    
penalty2 <- function(due, submit)
{
  due = strptime(due,  format = "%b %d %H:%M")
  sub = strptime(submit, format = "%b %d %H:%M")
  dayslate = as.numeric(difftime(sub, due, units = "days"))
  halflife = 7
  expshape = 1 
  round(exp( log(.5)/halflife^expshape*(dayslate)^expshape ), 2)
}

#=================================================================================================================================
                                    
pen.plot2 <- function(dayslate = 50){

submit <- strftime(seq(as.POSIXct("2020-07-23 10:20"), by = "1 day", length.out = dayslate+1), "%b %e %H:%M")
plot(0:dayslate, penalty("Jul 23 10:20", submit), type = "l", las = 1, tck = -0.03,
    xaxt = "n", xlab = "Days Late", ylab = "Penalty", lwd = 2, mgp = c(2, .4, 0), cex.axis = .9)
axis(1, at = 0:dayslate, cex.axis = .6,mgp = c(2, .01, 0), tck = -0.03)
}


#=================================================================================================================================
                                    
my.penalty2 <- function(dayslate = 0, dayslate.span = 30){
  
  dayslate.span <- tail(round(abs(dayslate.span)), 1)
  dayslate <- tail(round(abs(dayslate)), 1)
  if(dayslate > dayslate.span) dayslate.span <- dayslate
  pen.plot2(dayslate.span)
  submit <- strftime(seq(as.POSIXct("2020-07-23 10:20"), by = "1 day", length.out = dayslate+1), "%b %e %H:%M")
  if(dayslate != 0) submit <- submit[-seq_len(dayslate)]
  x <- dayslate
  y <- penalty2("Jul 23 10:20", submit)
  points(x, y, type = "h", col = if(dayslate != 0) 2 else 1)
  points(x, y, bg = 'cyan', col = 'magenta', pch = 21, cex = 1.5)
  text(x, y, y, cex = .75, font = 2, pos = 3, xpd = NA, col = if(dayslate != 0) 2 else 1)
} 
       
#=================================================================================================================================
       
penalty <- function(dayslate)
{
  halflife = 7
  expshape = 1 
  round(exp( log(.5)/halflife^expshape*(dayslate)^expshape ), 2)
}

#=================================================================================================================================

pen.plot <- function(dayslate = 50){
  
  curve(penalty(x), 0, dayslate, las = 1, tck = -0.03,
       xaxt = "n", xlab = "Days Late", ylab = "Penalty", lwd = 2, mgp = c(2, .4, 0), cex.axis = .9)
  axis(1, at = 0:dayslate, cex.axis = .6, mgp = c(2, .01, 0), tck = -0.03)
}

#=================================================================================================================================

my.penalty <- function(dayslate = 0, dayslate.span = 30){
  
  dayslate.span <- round(abs(dayslate.span))
  dayslate <- round(abs(dayslate))
  if(dayslate > dayslate.span) dayslate.span <- dayslate
  pen.plot(dayslate.span)
  x <- dayslate
  y <- penalty(dayslate)
  points(x, y, type = "h", col = ifelse(dayslate != 0, 2, 1))
  points(x, y, bg = 'cyan', col = 'magenta', pch = 21, cex = 1.5)
  text(x, y, y, cex = .75, font = 2, pos = 3, xpd = NA, col = ifelse(dayslate != 0, 2, 1))
}        

#=================================================================================================================================
       
sim_piece <- function(nPart = 100, G.r = .3, G.sds = c(2, 1, 2), e = .1,
                      betas = 2*c(10, 0, 0, 0, 3, .3, .6, .5, 2.7, .2, .3, 0, .7, 1),
                      type = c("normal", "logistic", "poisson"), seed = NULL){
  
  type <- match.arg(type) 
  set.seed(seed)
   
  # TRUE MODEL: y ~ (Time1+Time2) * group * engage + prof + (Time1 + Time2 | id) 
  ##Betas: 
  #[1] "(Intercept)"         "Time1"               "Time2"               "groupT"             
  #[5] "engage"              "profMed"             "profAdv"             "Time1:groupT"       
  #[9] "Time2:groupT"        "Time1:engage"        "Time2:engage"        "groupT:engage"      
  #[13] "Time1:groupT:engage" "Time2:groupT:engage"
  
  # C flat at Time1 goes down Time2, but T goes up by .5 and 
  # at time 2 it goes up by .7
  
  cor2cov <- function (sds, R) outer(sds, sds) * R
  
  Sigma <- cor2cov(G.sds, diag(1-G.r, 3)+G.r) # 2 time pieces and 1 id as random effects
  
  b <- MASS::mvrnorm(nPart, mu = rep(0, nrow(Sigma)), # random deviation from average model
                     Sigma = Sigma)                   # Equivalent of betas in fixed-effects
  
  nTime <- 6 # Number of total measurements per person
  
  DF <- data.frame(id = id <- rep(1:nPart, each = nTime),
                   Time  = rep(0:5, nPart),
                   Time1 = c(0, 1, 2, 2, 2, 2),
                   Time2 = c(0, 0, 0, 1, 2, 3),
                   group = factor(sample(c("C", "T"), nPart, replace = TRUE)[id], levels = c("C", "T")),
                    prof = factor(sample(c("Beg", "Med", "Adv"), nPart, replace = TRUE)[id], 
                                 levels = c("Beg", "Med", "Adv")),  # Putting levels takes first act. as ref.
                   engage = rbeta(nPart, 2, 5)[id]) 
  
  # fixed-effect design matrix
  X <- model.matrix(~ (Time1+Time2) * group * engage + prof, data = DF)
  
  # Random-effects design matrix
  Z <- model.matrix(~ Time1+Time2, data = DF) 
  
  # linear predictors (fixed + random):
  lin.pred <- as.vector(X %*% betas + rowSums(Z * b[id,]))
  
  # outcome y:
  DF$y <- switch(type, normal = rnorm(nPart * nTime, lin.pred, e), 
                 logistic = rbinom(nPart * nTime, 1, plogis(lin.pred)), 
                 poisson = rpois(nPart * nTime, exp(lin.pred)))
  
  return(DF)
}       

#=================================================================================================================================
       
pt.curve <- function(X, adjust = 1, compact = NULL, pch = 16, col = 2, cex = .7, seed = 0, reset = TRUE, add = FALSE, na.rm = TRUE, ...) {
  
  if(na.rm) X <- na.omit(X)  
  n.target <- length(X)
  
  d <- density(X, adjust = adjust, n = n.target)
  
  n <- if(!is.null(compact)) { 
    
    auc <- sum(d$y*median(diff(d$x)))/(diff(range(d$x))*max(d$y))
    
    compact*ceiling(n.target/auc)
    
  } else { n.target }
  
  set.seed(seed)
  pts <- data.frame(x = runif(n, min(d$x), max(d$x)), y = runif(n, 0, max(d$y)))
  
  pts <- pts[pts$y < approx(d$x, d$y, xout = pts$x)$y, ]
  
  if(nrow(pts) == 0) stop("Increase the size of sample 'X' OR use 'compact = NULL'.", call. = FALSE)
  
  pts <- pts[sample(seq_len(nrow(pts)), n, replace = TRUE), ]
  
  if(!add){
  
  if(reset) graphics.off()    
  plot(pts, pch = pch, col = col, cex = cex, ...)
    
  } else {
    
  points(pts, pch = pch, col = col, cex = cex, ...)
    
  }
}             
       
#=================================================================================================================================
       
is.balanced <- function(formula, data, na.action){

!is.list(replications(formula, data, na.action))
}       

#=================================================================================================================================

na <- function(data, na_pat = '^[[:punct:]]+$'){

unique(grep(na_pat, as.character(unlist(data)), value = TRUE))

}       


#=================================================================================================================================         

total_sigma <- function(fit){ 
  
vc <- VarCorr(fit)
  
if(inherits(fit, "lme")) sqrt(sum(as.numeric(vc[,"Variance"]), na.rm = TRUE)) else sqrt(sum(as.numeric(c(attr(vc[[1]], "stddev"), attr(vc, "sc")))^2, na.rm = TRUE))
}     
       
#=================================================================================================================================         
       
t.testb <- function(m1, m2, s1, s2, n1, n2 = NA, m0 = 0, var.equal = FALSE, sdif = NA, r = NA, digits = 6){
  
  if(var.equal & !is.na(n2))
    {
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  } else if(!var.equal & !is.na(n2))
    {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    df <- ((s1^2/n1 + s2^2/n2)^2)/((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))
  }else
  {
    se <- if(!is.na(sdif)) sdif/sqrt(n1) else sdif(sdpre = s1, sdpos = s2, r = r)/sqrt(n1)
    df <- n1 - 1
  }
  
  t <- (m2-m1-m0)/se
  
  a <- round(data.frame(mean.dif = m2-m1, std.error = se, t.value = t, p.value = 2*pt(-abs(t),df)), digits)
  a$paired <- if(is.na(n2)) TRUE else FALSE
  a    
}              

#=================================================================================================================================
       
G_matrix2 <- function(fit){
  vc <- VarCorr(fit)
  if(length(vc)==1)
  {
    y <- as.matrix(Matrix::bdiag(vc))
  } else
  {
    z <- do.call(rbind,vc)
    y <- as.matrix(Matrix::bdiag(vc))
    dimnames(y)[[1]] <- rownames(z)
    dimnames(y)[[2]] <- rownames(z)
  }
  return(y)
}  
       
#=================================================================================================================================  
       
       
G_matrix <- function(fit, digits = 8){
  
  vc <- VarCorr(fit)
  
  if(inherits(fit, c("lmerMod", "lmerModLmerTest", "lme4"))){
    
    out <- as.matrix(Matrix::bdiag(vc))
    if(is.null(unlist(dimnames(out)))) {
      nm <- unlist(lapply(vc, function(x) attributes(x)$dimnames[1]))
      dimnames(out) <- list(nm, nm)
    }
    round(out, digits)
    
  } else if(inherits(fit, "lme") & length(fit$group) < 2) { round(getVarCov(fit), digits) } else { vc }
  
}       
       
   
#=================================================================================================================================  
    
has_warning <- function(m) {
  df <- summary(m)
  !is.null(df$optinfo$conv$lme4$messages) && 
    grepl('failed to converge', df$optinfo$conv$lme4$messages)
}

#===============================================================================================================================

converge1 <- function(fit, parallel = c("multicore","snow","no")[3], maxfun = 1e5){
  
  if(has_warning(fit)){
    diff_optims <- lme4::allFit(fit, maxfun = maxfun, parallel = parallel)
    is.OK <- sapply(diff_optims, is, "merMod")
    diff_optims.OK <- diff_optims[is.OK]
    
    convergence_results <- lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)
    working_indices <- sapply(convergence_results, is.null)
    if(sum(working_indices)==0){
      print("No algorithms from allFit converged.")
      print("You may still be able to use the results, but proceed with extreme caution.")
      first_fit <- NULL
    } else {
      first_fit <- diff_optims[working_indices][[1]]
    }
    if(isSingular(first_fit, tol = 1e-3))  message("The model is singular.") else first_fit
    
  } else if(isSingular(fit, tol = 1e-3)){ message("The model has converged but is singular.") 
  } else { message("No issues found with the model.")}
}
                                
#===============================================================================================================================


converge2 <- function(fit){
  
  if(has_warning(fit)){
    optimx_options <- c("L-BFGS-B", "nlminb", "nlm", "bobyqa", "nmkb", "hjkb")
    
    for(i in 1:length(optimx_options)){
      model_flex <- suppressWarnings(update(fit,  
                           control = lmerControl(optimizer = "optimx", optCtrl = list(method = optimx_options[i]))))
      
      if(is.null(model_flex@optinfo$conv$lme4$messages)){
        print(paste0("The optimx option '", optimx_options[i],"' worked!"))
        return(model_flex)
        
        if(isSingular(model_flex, tol = 1e-3)) message("The model has converged but is singular.")
        
        break
      } else { print(paste0("The optimx option '", optimx_options[i],"' didn't work"))}
    }
  }
  
  else if(isSingular(fit, tol = 1e-3)){ message("The model has converged but is singular.") 
  } else { message("No issues found with the model.")}
  
}
#===============================================================================================================================

converge3 <- function(fit){
  
  if(has_warning(fit)){
    algoptions <- c("NLOPT_LN_PRAXIS", "NLOPT_GN_CRS2_LM",
                    "NLOPT_LN_COBYLA", "NLOPT_LN_NEWUOA",
                    "NLOPT_LN_NEWUOA_BOUND", "NLOPT_LN_NELDERMEAD",
                    "NLOPT_LN_SBPLX", "NLOPT_LN_BOBYQA")
    
    for(i in 1:length(algoptions)){
      
      model_flex <- suppressWarnings(update(fit, control = lmerControl(optimizer = "nloptwrap", optCtrl = list(method = algoptions[i],
                                                                                              maxit = 1e9,
                                                                                              maxeval = 1e9,
                                                                                              maxfun = 1e9,
                                                                                              xtol_abs = 1e-9,
                                                                                              ftol_abs = 1e-9))))
      if(is.null(model_flex@optinfo$conv$lme4$messages)){
        print(paste0("The nloptwrap options '", algoptions[i],"' worked!"))
        return(model_flex)
        break
      } else { print(paste0("The 'nloptwrap' options '", algoptions[i],"' didn't work!"))}
    }
    
  } else if(isSingular(fit)){ message("The model has converged but is singular.") 
  } else { message("No issues found with the model.")}    
  
}
#==================================================================================================================================


par_restart <- function(fit){
  
  if(has_warning(fit)){
    strict_tol <- lmerControl(optCtrl=list(xtol_abs =1e-8, ftol_abs=1e-8))  
    
    if (isLMM(fit)) {
      pars <- getME(fit,"theta")
    } else {
      pars <- getME(fit, c("theta","fixef"))
    }
    
    restart <- suppressWarnings(update(fit, start = pars))
    
    restart1 <- if(!is.null(restart@optinfo$conv$lme4$messages)) NULL else restart
    
    mins <- pmin(pars/1.01, pars*1.01)
    maxs <- pmax(pars/1.01, pars*1.01)
    
    pars_x <- runif(length(pars), mins, maxs) 
    
    restart2 <- suppressWarnings(update(fit, start=pars_x,
                       control=strict_tol))
    
    restart2 <- if(!is.null(restart2@optinfo$conv$lme4$messages)) NULL else restart2
    
    if(is.null(restart1)) message("Note: 'restart1' didn't work (hence NULL)!") 
    if(is.null(restart2)) message("Note: 'restart2' didn't work (hence NULL)!")
    
    list(restart1 = restart1, restart2 = restart2)
    
  } else if(isSingular(fit)){ message("The model has converged but is singular.") 
  } else { message("No issues found with the model.")}                                 
  
}
#=================================================================================================================================


par_compare <- function(fit){
  
  if(has_warning(fit)){
  fit.all <- lme4::allFit(fit)
  ss <- summary(fit.all)
  list(fixef = ss$ fixef,            ## fixed effects
       logLike = ss$ llik,                ## log-likelihoods
       SDs_Cor = ss$ sdcor,               ## SDs and correlations
       Cholesky = ss$ theta,               ## Cholesky factors
       did_algorithms_work_OK = ss$ which.OK)     ## which fits worked
  
  } else if(isSingular(fit)){ message("The model has converged but is singular.") 
  } else { message("No issues found with the model.")}                                

}

#=================================================================================================================================
                                
                                
par_hess_grad <- function(fit){

devfun <- update(fit, devFunOnly = TRUE)

if (isLMM(fit)) {
  pars <- getME(fit,"theta")
} else {
  pars <- getME(fit, c("theta","fixef"))
}
if (require("numDeriv")) {
  hessianA <- hessian(devfun, unlist(pars))
  gradientA <- grad(devfun, unlist(pars))
}
lme4_cal <- fit@optinfo$derivs

list(lme4_hessian = lme4_cal$Hessian,external_hessian = hessianA, lme4_gradient = lme4_cal$gradient, 
external_gradient = gradientA)
}                                
                                
                                                                                        
#=================================================================================================================================
                                
                                
find.norms <- function(low, high, cover = .99, digits = 6){
   
   f <- Vectorize(mu.norm)
   data.frame(t(f(low = low, high = high, cover = cover, digits = digits)))
 }              
                
#====================================================================================================
                
mu.norm <- find.norm <- function(low, high, cover = .99, digits = 6){
  
  options(warn = -1)
  
  cover[cover >= 1] <- .999999999999999
  cover[cover <= 0] <- .000000000000001
  
  p1 <- (1 - cover) / 2 
  p2 <- 1 - p1
  
  q <- c(low, high)  
  alpha <- c(p1, p2)
  
  is.df <- function(a, b, sig = 4) (round(a, sig) != round(b, sig))
  
  if (p1 <= 0 || p2 >= 1 || q[1] >= q[2] || p1 >= p2) {
    
    stop("Incorrect 'low' and/or 'high' or 'cover' values.", call. = FALSE)
    
  } else {
    
    beta <- qnorm(alpha)
    
    parm <- solve(cbind(1, beta), q)
    
    q <- qnorm(c(p1, p2), parm[[1]], parm[[2]])
  }
  
  if(is.df(low, q[[1]]) || is.df(high, q[[2]])) {
    
    stop("Change 'low' and/or 'high' or 'cover' values.", call. = FALSE)
    
  } else {
    
    return(round(c(mean = parm[[1]], sd = parm[[2]]), digits = digits))
  }
}                                                 
                                
#========================================================================================================================================================================
                                
                                
cor2G <- function(sd.int = 6, sd.slope = .01, rho = .3){
  
cormat <-  matrix(c(sd.int, rho, rho, sd.slope), 2, 2)
res <- data.frame(lme4::sdcor2cov(cormat), row.names = c("Int.", "slope") ) 
colnames(res) <- rownames(res)
as.matrix(res)
}


#==========================================================================================================================================================================================


sim_m4b <- function(n_cluster, ave_cluster_n, G = cor2G(sd.int = 2, sd.slope = .275, rho = 1), G.r = .3, G.sds = c(6, .01), betas = c(12, 3, 2, -1.3),
                         e = 3.2, empirical = TRUE, seed = NULL, output_data = FALSE, re.var = "sd"){

# math ~ ses * sector + (ses | sch.id)  
set.seed(seed)

auto <- missing(n_cluster) | missing(ave_cluster_n)  
  
if(auto){
hsb <- read.csv('https://raw.githubusercontent.com/rnorouzian/e/master/hsb.csv')
n_cluster <- 160
ave_cluster_n <- hsb %>% count(sch.id) %>% pull(n)
}
  
cor2cov <- function (sds, R) outer(sds, sds) * R

#G <- cor2cov(G.sds, diag(1-G.r, length(G.sds))+G.r)

b <- MASS::mvrnorm(n_cluster, mu = rep(0, nrow(G)), 
                   Sigma = G, empirical = empirical)                   

#sector = factor(hsb$sector)
#factor(sample(c("pub", "cath"), n_cluster, replace = TRUE)[sch.id], levels = c("pub", "cath"))

DF <- data.frame(sch.id = sch.id <- rep(1:n_cluster, if(auto) ave_cluster_n else each = ave_cluster_n), 
                  sector = factor(hsb$sector),
                  ses = rnorm(n_cluster, -.5, .65)[sch.id])

data_size <- nrow(DF)

# fixed-effect design matrix
X <- model.matrix(~ ses*sector, data = DF)

# Random-effects design matrix
Z <- model.matrix(~ ses, data = DF)

# linear predictors (fixed + random):
lin.pred <- as.vector(X %*% betas + rowSums(Z * b[sch.id,]))

# outcome math:
DF$math <- rnorm(data_size, lin.pred, e)

fit <- lmer(math ~ ses * sector + (ses | sch.id), data = DF)
G <- G_matrix(fit)

list(fit = fit, G_matrix = G, G_cor_matrix = stats::cov2cor(G), PCA = G_pca(fit), data = if(output_data) DF else NULL)

}
          
     
#===================================================================================================================================     
     
sim_m4 <- function(n_cluster, ave_cluster_n, G = cor2G(sd.int = 2, sd.slope = .275, rho = 1), G.r = .3, G.sds = c(6, .01), betas = c(12, 3, 2, -1.3),
                         e = 3.2, empirical = TRUE, seed = NULL, output_data = FALSE, re.var = "sd"){

# math ~ ses * sector + (ses | sch.id)  
set.seed(seed)

auto <- missing(n_cluster) | missing(ave_cluster_n)  
  
if(auto){
hsb <- read.csv('https://raw.githubusercontent.com/rnorouzian/e/master/hsb.csv')
n_cluster <- 160
ave_cluster_n <- hsb %>% count(sch.id) %>% pull(n)
}
  
cor2cov <- function (sds, R) outer(sds, sds) * R

#G <- cor2cov(G.sds, diag(1-G.r, length(G.sds))+G.r)

b <- MASS::mvrnorm(n_cluster, mu = rep(0, nrow(G)), 
                   Sigma = G, empirical = empirical)                   

#sector = factor(hsb$sector)
#sector = factor(sample(c("pub", "cath"), n_cluster, replace = TRUE)[sch.id], levels = c("pub", "cath"))

DF <- data.frame(sch.id = sch.id <- rep(1:n_cluster, if(auto) ave_cluster_n else each = ave_cluster_n), 
                  sector = factor(hsb$sector),
                  ses = rnorm(if(auto) nrow(hsb) else n_cluster*ave_cluster_n, -.5, .65))

data_size <- nrow(DF)

# fixed-effect design matrix
X <- model.matrix(~ ses*sector, data = DF)

# Random-effects design matrix
Z <- model.matrix(~ ses, data = DF)

# linear predictors (fixed + random):
lin.pred <- as.vector(X %*% betas + rowSums(Z * b[sch.id,]))

# outcome math:
DF$math <- rnorm(data_size, lin.pred, e)

fit <- lmer(math ~ ses * sector + (ses | sch.id), data = DF)

return(fit)

}     
     
#=================================================================================================================================
     

G_pca <- function(fit) {
  
  obj <- summary(lme4::rePCA(fit))
  model <- lme4::VarCorr(fit)
  if(length(obj) == length(model)) {
    obj <- Map(function(x, z) {
      colnames(x$importance) <- paste(z, colnames(model[[z]]), sep = '_')
      x
    }, obj, names(obj))
  }
  else if(length(obj) == 1) {
    colnames(obj[[1]]$importance) <- unlist(mapply(paste, names(model), sapply(model, colnames), MoreArgs = list(sep = '_')))
  }
  return(obj)
}

     
#=================================================================================================================================
     
scale_fit <- function(fit){
  
summary(jtools::scale_mod(fit))
  
}       

#=================================================================================================================================                              
                              
group.center <- function (var, grp) 
{
  grp <- as.factor(grp)
  grp <- as.numeric(grp)
  var <- as.numeric(var)
  return(var - tapply(var, grp, mean, na.rm = TRUE)[grp])
}

#===============================================================================================================================
                                      
group.mean <- function (var, grp) 
{
  grp <- as.factor(grp)
  grp <- as.numeric(grp)
  var <- as.numeric(var)
  return(tapply(var, grp, mean, na.rm = TRUE)[grp])
}   
                              
#=================================================================================================================================  
                 
                 
sim_sng <- function(n_cluster = 10, ave_cluster_n = 15, G.r = 0, G.sds = c(2, 1, 1, 1), betas = c(-10, 10, 10, 10),
                   e = 2, empirical = FALSE, seed = NULL, output_data = FALSE){
  
  # y ~ A + B + C + (A + B + C | group)  
  set.seed(seed)
  
  data_size <- n_cluster*ave_cluster_n
  
  cor2cov <- function (sds, R) outer(sds, sds) * R
  
  G <- cor2cov(G.sds, diag(1-G.r, length(G.sds))+G.r)
  
  b <- MASS::mvrnorm(n_cluster, mu = rep(0, nrow(G)), 
                     Sigma = G, empirical = empirical)                   
  
  DF <- data.frame(group = group <- rep(1:n_cluster, each = ave_cluster_n), 
                   A = rnorm(data_size, 10),
                   B = rnorm(data_size, 10),
                   C = rnorm(data_size, 10))
  
  # fixed-effect design matrix
  X <- model.matrix(~ A+B+C, data = DF)
  
  # Random-effects design matrix
  Z <- model.matrix(~ A+B+C, data = DF)
  
  # mu: linear predictors (fixed + random):
  mu <- as.vector(X %*% betas + rowSums(Z * b[group,]))
  
  # outcome math:
  DF$y <- rnorm(data_size, mu, e)
  
  return(round(DF, 2))
}  

#=================================================================================================================================
                              
R2_null <- function(m0, ..., digits = 4){

foo <- function(m0, m1){
  
  VarCorr(m0) %>% 
    as.data.frame %>% 
    select(grp, var_m0 = vcov) %>% 
    left_join(VarCorr(m1) %>% 
                as.data.frame %>% 
                select(grp, var_m1 = vcov), by = "grp") %>% 
    mutate(var_red = 1 - var_m1 / var_m0) 
}

temp <- lapply(list(...) , foo, m0 = m0)
temp2 <- cbind(lapply(temp, '[', 4))

result <- round(do.call(cbind, temp2), digits)
names(result) <-  substitute(...())
rownames(result) <- c("Level-2:", "Level-1:")
result
}                              

#=================================================================================================================================  
     
                              
do_context <- function(data, context_vars, group_id){
  
  all_names <- names(data)
  
  id <- grep("id|group|grp", all_names, value = TRUE, ignore.case = TRUE)
  
  if(!all(group_id %in% all_names)) { 
    
    stop(paste(toString(dQuote(group_id)), "not found for 'group_id' in the 'data'.", if(length(id)>0) 
      paste("\nPossibly you meant to use", toString(dQuote(id)), "as 'group_id', no?")), call. = FALSE) 
    
    }
  
  ok <- context_vars %in% all_names
  
  if(!all(ok)) message(paste("\n", toString(dQuote(context_vars[!ok])), "not found in the 'data' thus ignored."))
  
  context_vars <- context_vars[ok] 
  
  dum_vars <- all_names[sapply(data, function(i)is.character(i)|is.factor(i))]
  
  dum_names <- context_vars[context_vars %in% dum_vars]
  
  is_dum <- length(dum_names) > 0
  
  num_names <- context_vars[!(context_vars %in% dum_vars)]
  
  is_num <- length(num_names) > 0
  
  
  if(is_num){
    data <- data %>%
      group_by(across(all_of(group_id))) %>%                          
      mutate(across(all_of(num_names), list(wthn = ~ . - mean(.), btw = ~ mean(.)))) %>% 
      as.data.frame()
  }
  
  
  if(is_dum){
    data <- data %>%
      dummy_cols(select_columns = dum_names) %>% 
      group_by(across(all_of(group_id))) %>%                          
      mutate(across(starts_with(paste0(dum_names, "_")), list(wthn = ~ . - mean(.), btw = ~ mean(.)))) %>% 
      as.data.frame()
  }
  
  return(data)
}

#=================================================================================================================================  
                               
                               
hetro_var <- function(fit){
  
    if(!inherits(fit, c("lme", "gls"))) stop("Only 'lme()' & 'gls()' models are accepted.", call. = FALSE)
  
  coef(fit$modelStruct$varStruct, uncons = FALSE, allCoef = TRUE)
  
}
     
#=================================================================================================================================

                               
rho_lme <- function(fit) {
  
  if(!inherits(fit, c("lme", "gls"))) stop("Only 'lme()' & 'gls()' models are accepted.", call. = FALSE)
  
  coef(fit$modelStruct$corStruct, uncons = FALSE, allCoef = TRUE)

}                               
 
                               
#=================================================================================================================================
                               
cov_str <- function(fit, cov = TRUE, time_var = "time", hlm = TRUE){
  
  rho <- rho_lme(fit)
  hetro <- hetro_var(fit)
  sig <- sigma(fit)
  dat <- getData(fit)
  if(!(time_var %in% names(dat))) stop("Your 'time_var' doesn't exist in your data.", call. = FALSE)
  time_vals <- unique(dat[[time_var]])
  
  if(is.null(rho) & is.null(hetro)) return(id_cor(fit, cov = cov, time_var = time_var, hlm = hlm))
  if(is.null(rho) & !is.null(hetro)) {
    
    res <- id_cor(fit, cov = cov, time_var = time_var, hlm = hlm)
    diag(res) <- sum_ranef_var(fit)+(sig*hetro)^2
    return(res)  
    
  }
  corm <- corMatrix(fit$modelStruct$corStruct)[[1]]
  
  res <- corm*sig^2+if(hlm)sum_ranef_var(fit) else 0 *if(is.null(hetro)) 1 else t(t(hetro))%*%t(hetro)
  if(!is.null(hetro)) diag(res) <- sum_ranef_var(fit)+(sigma(fit)*hetro_var(fit))^2
    
  if(!cov) res <- cov2cor(res)
  rownames(res) <- colnames(res) <- paste0(time_var,time_vals)
  return(res)
  
  }                               
                               
#=================================================================================================================================  
                               
id_cor <- function(fit, cov = TRUE, time_var = "time", hlm = TRUE){
  
  sig <- sigma(fit)^2
  vrs <- as.numeric(VarCorr(fit)[,"Variance"])
  time_vals <- unique(getData(fit)[[time_var]])
  steps <- length(time_vals)
  x <- diag(steps)
  res <- sig*x
  
  diag(res) <- sum_ranef_var(fit, resid = TRUE) 
  res[col(res)!=row(res)] <- if(hlm) sum_ranef_var(fit) else 0
  
  rownames(res) <- colnames(res) <- paste0(time_var,time_vals)
  if(!cov) res <- cov2cor(res)
  return(res)
} 
                               
#=================================================================================================================================
           
sum_ranef_var <- function(fit, resid = FALSE){
vrs <- as.numeric(VarCorr(fit)[,"Variance"])
sum(if(resid) vrs else rev(vrs)[-1], na.rm = TRUE)
}                               
                               
#=================================================================================================================================                               
                               
quad <- function(x, a, b, c, degree = 2) a + b*x + c*x^degree   
                               
#=================================================================================================================================                               
                               
get_forms <- function(dredge_fit, n = 1:5){ 
  
  zz <- get.models(dredge_fit, n)
  
  data.frame(model = as.character(lapply(seq_along(zz), 
                                         function(i)get_formula(zz[[i]]))))
  }                               
                               
#=================================================================================================================================  
  
need <- c("lme4", "nlme", "glmmTMB", "emmeans", "plotrix", "ellipse", 'jtools', 'stargazer', 'interactions', 'car', 'MASS', 'modelr', 'fastDummies', 'MuMIn', 'sjPlot', 'lmerTest', 'reghelper',
          'bbmle', 'performance', 'see', 'psych','haven', 'effects','tidyverse','parallel','optimx','minqa','blme','dfoptim', 'remotes', 'DHARMa', 'multcomp', 'splines')
     
not.have <- need[!(need %in% installed.packages()[,"Package"])]
if(length(not.have)) install.packages(not.have)

#options(warn = -1)

suppressWarnings(                                         
suppressMessages({ 
  
for(i in need){
  library(i, character.only = TRUE)
}
}))                                
