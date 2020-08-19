
Break = "\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = "    Welcome to 'EDP 380C.16 Hierarchical Linear Modeling'.
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

plot.cor <- function (corr, outline = FALSE, col = colorRampPalette(c(4, 2))(choose(ncol(corr), 2)), 
                      upper.panel = c("ellipse", "number", "none"), lower.panel = c("ellipse", "number", "none"), diag = c("none", "ellipse", "number"), digits = 2, bty = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, cex.lab = par("cex.lab"), cex = 0.75 * par("cex"), mar = c(2, 2, 4, 2)-.2, ...)
{
  savepar <- par(pty = "s", mar = mar)
  on.exit(par(savepar))
  if (is.null(corr))
    return(invisible())
  if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1))
    stop("Please input a correlation matrix")
  plot.new()
  par(new = TRUE)
  rowdim <- dim(corr)[1]
  coldim <- dim(corr)[2]
  rowlabs <- dimnames(corr)[[1]]
  collabs <- dimnames(corr)[[2]]
  if (is.null(rowlabs))
    rowlabs <- 1:rowdim
  if (is.null(collabs))
    collabs <- 1:coldim
  rowlabs <- as.character(rowlabs)
  collabs <- as.character(collabs)
  col <- rep(col, length = length(corr))
  dim(col) <- dim(corr)
  upper.panel <- match.arg(upper.panel)
  lower.panel <- match.arg(lower.panel)
  diag <- match.arg(diag)
  cols <- 1:coldim
  rows <- 1:rowdim
  maxdim <- max(length(rows), length(cols))
  plt <- par("plt")
  xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", cex = cex.lab))/(plt[2] - plt[1])
  xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
  ylabwidth <- max(strwidth(collabs[cols], units = "figure", cex = cex.lab))/(plt[4] - plt[3])
  ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
  plot(c(-xlabwidth - 0.5, maxdim + 0.5), c(0.5, maxdim + 1 + ylabwidth), type = "n", bty = bty, axes = axes, xlab = "", ylab = "", asp = asp, cex.lab = cex.lab, ...)
  text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], adj = 1, cex = cex.lab)
  text(cols, rep(length(rows) + 1, length(cols)), labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
  mtext(xlab, 1, 0)
  mtext(ylab, 2, 0)
  mat <- diag(c(1, 1))
  plotcorrInternal <- function() {
    if (i == j){ 
      if (diag == 'none'){
        return()
      } else if (diag == 'number'){
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else if (diag == 'ellipse') {
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      }
    } else if (i >= j){ 
      if (lower.panel == 'ellipse') { 
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (lower.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    } else { 
      if (upper.panel == 'ellipse') { 
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (upper.panel == 'number') { 
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    }
  }
  for (i in 1:dim(corr)[1]) {
    for (j in 1:dim(corr)[2]) {
      plotcorrInternal()
    }
  }
  invisible()
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
                                                                       
lo_ave_up <- function(data, vars) sapply(vars, function(x) 
  setNames(mean(data[[x]]) + c(-1, 0, 1)*sd(data[[x]]), 
           paste0(x, c('-1SD', '.Mean', '+1SD'))), simplify = FALSE)                                    
                                                                        
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
  
  R <- attr(lme4::VarCorr(fit)[[group.id]], "correlation")
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
    SDs <- diag(attr(lme4::VarCorr(fit)[[group.id]], "stddev"))
    (C <- SDs %*% R %*% SDs)
    res <- base::solve(C[pred_loc, pred_loc], C[out_loc, pred_loc], tol = tol)
  }
  names(res) <- if(length(pred.names) > 1) colnames(R[pred_loc, pred_loc]) else pred.names
  return(res)
}


#=================================================================================================================================================


latent.lmer <- function (fit, formula, group.id, std = TRUE, digits = 3, prog.bar = "none",
                             sing.check = TRUE, ridge.constant = 1e-06, tol = 1e-12,
                             seed = 123, nsim = 499, level = .95, parallel = c("no", "multicore", "snow")[2], ncpus = 7) {
  
  formals(latent.reg)[formalArgs(latent.reg)[-1L]] <- c(
    formula, group.id, std, sing.check, ridge.constant, tol)
  
  boot.res <- lme4::bootMer(
    fit, nsim = nsim, seed = seed, .progress = prog.bar,
    parallel = parallel, ncpus = ncpus,
    FUN = latent.reg)
  
  #if(inherits(boot.res, "try-error")) stop("Change 'parallel' and/or 'ncpus' values.")
  
  se.s <- apply(boot.res$t, 2, sd, na.rm = TRUE)
  z.value <- abs((2 * boot.res$t0 - colMeans(boot.res$t, na.rm = TRUE)) / se.s)
  
  boot.p <- 2 * pnorm(-z.value)
  
  ci <- confint(boot.res, level = level, type = "norm")
  
  round(data.frame(est. = boot.res$t0, SE = se.s, z.value = z.value, p.value = boot.p,
    lower = ci[[1]], upper = ci[[2]]), digits)
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
  random <- length(getME(fit,"theta"))
  error <- if(isLMM(fit)) 1 else 0
  return(c(total = total, fixed = fixed, random = random, error = error))
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
       
sim.piece <- function(nPart = 100, G.r = .3, G.sds = c(2, 1, 2), e = .1,
                      betas = 2*c(10, 0, 0, 0, 3, .3, .6, .5, 2.7, .2, .3, 0, .7, 1),
                      type = c("normal", "logistic", "poisson"), seed = NULL){
  
  type <- match.arg(type) 
  set.seed(seed)
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
  
need <- c("lme4", "nlme", "glmmTMB", "emmeans", "plotrix", "ellipse", 'jtools', 'stargazer', 'interactions', 'car', 'tidyverse', 'modelr', 'bbmle', 'performance', 'see','MASS', 'psych','haven', 'effects') 
not.have <- need[!(need %in% installed.packages()[,"Package"])]
if(length(not.have)) install.packages(not.have)

options(warn = -1)
       
suppressMessages({ 
  
for(i in need){
  library(i, character.only = TRUE)
}
})                                
