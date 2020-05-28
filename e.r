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

post.mixed <- function(fit, formula = NULL, plot = TRUE, by = NULL, var = NULL, horiz = TRUE, adjust = "tukey", type = "response", compare = FALSE, ...){
  
  limit <- nobs(fit)
  tm <- terms(fit)
  
  f <- if(is.null(formula)) as.formula(bquote(pairwise ~ .(tm[[3L]]))) else as.formula(formula)
  
  av <- emmeans::.all.vars(tm)
  
  cl <- if(inherits(fit, "lme")) attr(tm, "dataClasses")[-1L] else setNames(sapply(av[-1L], function(i) class(model.frame(fit)[[i]])), av[-1L])
  
  all.factor <- all(cl == "factor") || all(cl == "character") || all(cl == "character" | cl == "factor")
  
  no.contrast <- length(av[-1L]) == 1L & !all.factor
  
  ems <- if(all.factor || no.contrast)  { eval(substitute(emmeans::emmeans(fit, f, infer = c(TRUE, TRUE), type = type, pbkrtest.limit = limit))) 
    
  } else {
    
    var <- if(is.null(var)) names(sort(cl[grep("integer|numeric", cl)])[1]) else var
    
    f <- rm.terms(f, var)
    
    eval(substitute(emmeans::emtrends(fit, f, var = var, infer = c(TRUE, TRUE), type = type, pbkrtest.limit = limit)))
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
  names(res)[c(2,5:7, 9:11)] <- c(if(all.factor)"mean.dif"else paste0("slope.dif","(", var,")"), "lower", "upper", "t.value", "Cohen.d", "lower.d", "upper.d")
  res
  
  }
  
  else {
    
    ems
    
  }
  
  return(out)
}      
                                  
#=================================================================================================================================
                                  
data.tree <- function(data, toplab = NULL, cex = 1, auto = FALSE, ...){
  
  if(auto){    
    cats <- sapply(data, Negate(is.numeric))  
    data <- data[cats]
    
    what <- names(data)
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
  
#=================================================================================================================================  
  
need <- c("lme4", "nlme", "glmmTMB", "emmeans", "plotrix", "ellipse", "ggplot2") 
have <- need %in% rownames(installed.packages())
if(any(!have)){ install.packages( need[!have] ) }

options(warn = -1)
suppressMessages({ 

  library('lme4')
  library('nlme')
  library('glmmTMB')
  library('emmeans')
  library('ggplot2')
  library('plotrix')
  library('ellipse')
})  
