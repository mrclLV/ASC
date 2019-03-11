library(MASS)
library(Matrix)

Rho <- function(dim, uptri){
  S <- diag(dim)
  S[lower.tri(S)] <- uptri
  S <- S + t(S)
  S <- S - diag(dim)
  S
}


approxqf <- function(x){
  cdf <- ecdf(c(x))
  smoothcdf <- approxfun(x = x, cdf(x), yleft = 0, yright = 1)

  mnx <- min(x)
  mxx <- max(x)

  qd <- function(q){
    if(q <= 0.01) return(mnx)
    if(q >= 1) return(mxx)
    uniroot(function(x){ smoothcdf(x) - q}, interval = c(mnx, mxx) )$root
  }

  Vectorize(qd)

}

runifcor <- function(nms, n, rho){
  rho <- nearPD(rho)$mat
  lnms <- length(nms)

  Xnorm <- mvrnorm(n,mu =rep(0,lnms), rho)
  Xunif <- pnorm(Xnorm)


  X <- data.frame(Xunif)
  X
}




factor_runifcor <- function(x, lev, p = NULL){

  if(is.null(p)){

    brk <- lapply(lev, function(x){
                    if(any(x > 0)){
                      seq(0,1, by = 1/length(x))
                  } else {
                    0
                  }
                  }   )

  } else {
    brk <- lapply(p, function(x) c(0,cumsum(x)))
  }

  res <- Map( function(x, b) {
               if( any(b > 0)  ){

                 cut(x= x, breaks = b)

               } else {

                 x

               }},

             x = x, b= brk)

  ## res <- lapply(res, function(x){
  ##          factor(x, sample(levels(x)))
  ##        }    ) # shuffle to "avoid" linear behaviour


  res <- Map(`levels<-`, res, lev)

  data.frame(res)

}



cor_sim <- function(adjs, rho ,p = NULL, n = NULL){


  covrts <- Covariates(adjs)
  c_covrts <- Covariates(continous(adjs))
  whichc <- whichcont(adjs)


  lev <- levels(adjs)
  nms <- names(lev)

  if( is.null(n) ) n <- nrow(covrts)

  cor_unif <- runifcor(nms, n , rho)

  if( is.null(p) ){
    p <- lapply(covrts, table)
    p <- lapply(p,prop.table)

  } else {
    ep <- lapply(covrts, table)
    ep <- lapply(ep,prop.table)

    p_null <- vapply(p, is.null, logical(1))
    p[p_null] <- ep[p_null]
    p <- lapply(p, as.numeric)
  }


  cat_data <- factor_runifcor(cor_unif, lev, p)


  if(length(whichc) > 0 ){

    qfuncs <- lapply(c_covrts[whichc],approxqf)
    cont_unif <- lapply(cor_unif[whichc], unlist)
    cont_unif <- unname(cont_unif)
    cont_data <- Map(do.call, what =qfuncs,
                     args = list(cont_unif) )
     cont_data <- data.frame(cont_data)
  }

  names(cat_data) <- attr(terms(adjs),'term.labels')
  data.frame(cat_data,cont_data)
}



