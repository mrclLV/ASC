# Re - implementation of adjusted survival curves, with a more polished object oriented structure, so its easier to understand and more short to write. It tries to mimic the operation of survival package, with an AdjSurv as the central  object.  Also an important point, numerical variables may be used where possible (logistic regression, random forest, cox models..) instead of grouped variables. Error return control.

# And althought still functional with time_grid, by default, all comparisons wil be done at points where we have information.

library(survival)
library(randomForest)

## Auxiliar functions to use survival objects



scurve_survfit <- function(sfit, time){
  strata <- names(sfit$strata)
  strata <- gsub("(( ?, +)|( +, ?))",
                 strata ,replacement = '%&')
  strata <- gsub("\\w+=",
                  strata ,replacement = '')
  strata <- gsub(" ",
                 strata ,replacement = '')

  names(sfit$strata) <- strata
  sm <- summary(sfit, times  = time, extend = TRUE)

  strata <- as.character(sm$strata)
  groups <- strsplit(x =strata, split = '%&')

  groups <- do.call(rbind,groups)

  ## start <- data.frame(time = 0,
  ##                     surv = 1,
  ##                     se = 0,
  ##                     group =unique(unlist(groups))
  ##                     )

 # names(start)[4] <- 'group'
  ret <- data.frame(time = sm$time, surv = sm$surv, se = sm$std.err, group = groups)
  names(ret)[4] <- 'group'
  ret
}



scurve_survexp <- function(sexp){
  covariate_name <- terms(formula(sexp$call))
  covariate_name <- attributes(covariate_name)$term.labels
  lev <- colnames(sexp$surv)
  lev <- gsub(paste0(covariate_name,"="),lev ,replacement = '')
  group <- rep(lev, each = length(sexp$time))
  names(group) <- 'group'

  if (any(sexp$surv ==  0) ) {
    sexp$surv <- apply(sexp$surv, 2, function (x) {
                     # Avoid expected survival starts with
                         i <- 1    # 0, when no individuals n a group
                         while(x[i]  == 0 ) {
                           x[i] <- 1
                           i <- i +1
                         }
                         x
                       }
                       )
  }



  ## start <- data.frame(time = 0,
  ##                     surv = 1,
  ##                     group =unique(group)
  ##                     )

  ret <- data.frame(time = sexp$time,
                    surv = c(sexp$surv),
                    group,
                    row.names = NULL)

   ret
}


Sfit <- function(survobj, intrst,  weights =NULL){
  survfit(survobj ~ intrst, weights = weights)
}


.limit_weights <- function(w, weight_lim){

  if(is.numeric(weight_lim) &&
     length(weight_lim)  == 2 &&
     weight_lim[2] > weight_lim[1] ){

    weights <- vapply(w, function(x) min(x, weight_lim[2]),
                      numeric(1))
    weights <- vapply(weights, function(x) max(x, weight_lim[1]),
                      numeric(1))
    weights

  } else {
    warning('Use proper limits, c(lower bound, upper bound).')
  }

}

###

Sobj <- function(x){
  UseMethod("Sobj")
}

## AdjSurv object functions

cont <- function(cont, grouped){

    attr(grouped,'cont') <-  cont
    grouped

  }

interest <- function(x){ x}

AdjSurv <- function(formula, survobj,  data){

  tr <- terms(formula, specials = ('cont'))
  cont_idx <- attr(tr,'specials')$cont

  if(is.null(cont_idx)){

    attr(tr,'term.labels') <- all.vars(formula)
    mdl <- model.frame(tr, data = data)

    names(mdl) <- all.vars(formula)
    ret <- list(covariates = mdl, Surv = survobj)
    class(ret) <- c(class(ret), "AdjSurv")

  } else {

    attr(tr,'term.labels') <- all.vars(formula)[-cont_idx]
    mdl <- model.frame(tr, data = data)
    cont_names <- all.vars(formula)[cont_idx]
    names(mdl) <- all.vars(formula)[- cont_idx]
    names(mdl)[cont_idx] <- cont_names
    ret <- list(covariates = mdl, Surv = survobj)
    class(ret) <- c(class(ret), "AdjSurv")

  }


  ret
}



print.AdjSurv <- function(x){
  print(x$covariates)
}

Interest <- function(adjs){
adjs$covariates[ ,  attr(terms(adjs$covariates), 'response')]
}

Adjvars <- function(adjs){
adjs$covariates[ , - attr(terms(adjs$covariates), 'response'), drop = FALSE]
}


Covariates <- function(adjs){
  ret <- adjs$covariates
  attr(ret, 'terms') <- NULL
  ret
}

time.AdjSurv <- function(adjs){
  adjs$Surv[ ,1]
}


Sobj.AdjSurv <- function(adjs){
  adjs$Surv
}

formula.AdjSurv <- function(adjs){
  formula(terms(adjs$covariates))
#  formula(adjs$covariates)
}

terms.AdjSurv <- function(adjs){
 terms(adjs$covariates)
}

whichcont <- function(adjs){
  attr(terms(adjs),'specials')$cont
}

levels.AdjSurv <- function(adjs){
  lapply(adjs$covariates, levels)
}



continous <- function(adjs){

  continous <- whichcont(adjs)

  if( !is.null(continous)  ){
    adjs$covariates[ ,continous] <-  lapply(continous,
                                function(y){
                                  attr(adjs$covariates[ , y], 'cont')
                                }  )
  }
  adjs
}

## adjsfit object functions


surv <- function(adjsfit){
  adjsfit$surv
}


time.adjsfit <- function(adjsfit){
  adjsfit$time
}


group <- function(adjsfit){
  adjsfit$group
}


se <- function(adjsfit){
  adjsfit$se
}

# adjsfit inicialization

.adjsfit <- function(s){
  if(inherits(s, 'try-error')) {
    s <- list()
    s$surv <- NA_integer_
    class(s) <- c('try-error', class(s))
  }
  class(s) <- c('adjsfit', class(s))
  s
}

.weighted_adjsfit <- function(weights, adjs, weight_lim, type, time = NULL){

  if(inherits(weights, 'try-error')){
    s <- NA_integer_
    return(.adjsfit(s))
  }

  if(is.numeric(weight_lim)){
    weights <- .limit_weights(weights, weight_lim)
  } else {

    weights
  }

   switch(type,

         weights = weights,

         surv = {
    intrst <- Interest(adjs)
    s <- Sfit(adjs$Surv, intrst, weights)
    if( is.null(time)) time <-  adjs$Surv[ ,1]
    s <- scurve_survfit(s, time)
    .adjsfit(s)
  }  )
}


# adjsfit output
##print.adjsfit <- function(adjsfit, ...){

## }

plot.adjsfit <- function(adjsfit,...){
  mtime <- max(adjsfit$time)
  plot(1, type = 'n',xlim = c(0,mtime), ylim = c(0,1), ...)
  adjsfit_group <- split(adjsfit, adjsfit$group)

  invisible(lapply(seq_along(adjsfit_group), function(i)lines( x = adjsfit_group[[i]]$time, y = adjsfit_group[[i]]$surv, type = 's', lty = i)))

}


lines.adjsfit <- function(adjsfit,...){
  adjsfit_group <- split(adjsfit, adjsfit$group)
  invisible(lapply(adjsfit_group, function(x)lines( x = x$time, y = x$surv, type = 's', ...)))
}


## Adjusting Methods



kaplan_meier <- function(adjs, time = NULL, ...){
  s <- try({
      intrst <- Interest(adjs)
      fit <- Sfit(adjs$Surv, intrst)
      if( is.null(time)) time <-  adjs$Surv[ ,1]
      s <- scurve_survfit(fit, time)
      s
    })
      .adjsfit(s)
}



## Marginal Methods

# Adjustment by reference population

reference_pop <- function(adjs, ref = NULL, type ='surv', time = NULL, weight_lim = NULL, ...){

  w <- try({
  if(is.null(ref)){
    ref <- table(Adjvars(adjs))
  }

  prob_covrts <- prop.table(table(Covariates(adjs)))
  prob_ref <- c(prop.table(ref))
  weights <- prob_ref/c(prob_covrts)

  groups <- interaction(Covariates(adjs))
  group_weight <- data.frame(group = levels(groups),
                             weight = weights)
  weights <- group_weight$weight[match(groups
                                      ,group_weight$group)]
      weights
    })
 .weighted_adjsfit(w , adjs, weight_lim, type, time)
}


# Adjustment by inverse probability weighting overall adjusting vars p


ipw <- function(adjs, ref = NULL, type ='surv', time = NULL, weight_lim = NULL, ...){

  w <- try({

      covrts <- Covariates(adjs)
      covrts_groups <- interaction(covrts)
      cond_p <- prop.table(table(covrts), length(covrts))
      weights <- 1/expand.grid(cond_p)[covrts_groups,]
      weights

    })

  .weighted_adjsfit(w, adjs, weight_lim, type, time)
}

# reference should be a dataframe with correct amount of each covariate

ipw_function <- function(ref){

  ref <- as.data.frame(prop.table(table(ref),4))

  function(adjs, type ='surv', time = NULL, weight_lim = NULL, ...){
    w <- try({

        covrts <- Covariates(adjs)
        lr <- length(ref)
        ref_p <- ref[ ,lr]

        ref_int <- interaction(ref[-lr])
        covrts_int <- interaction(covrts)

        weights <- unname( 1/ref_p[match(covrts_int, ref_int)])

        weights

    })

  .weighted_adjsfit(w, adjs, weight_lim, type, time)
}

}

# Adjustment by inverse probability weighting by logistic regression


ipw_logi <- function( adjs,formula = NULL, type ='surv', use_cont = TRUE, time = NULL, weight_lim = NULL, ...){

  w <- try({
  intrst <- Interest(adjs)

  if (use_cont){
    adjvars <- Adjvars(continous(adjs))

  } else {

    adjvars <- Adjvars(adjs)

  }

  if(is.null(formula)){
    nms <- names(adjvars)
    formula <- as.formula(paste0('1 ~ ', paste0(nms, collapse = ' + ')))
  }

  len <- length(intrst)
  resp_level <- vapply(levels(intrst),
                       function(x) x == intrst,
                       logical(len))


  model_level <- lapply(data.frame(resp_level),
                        function(x){
                          formula <- update.formula(formula, x  ~ . )
                          glm(formula, data= data.frame(x, adjvars) )

                                    }  )


  weights <- as.data.frame(lapply(model_level, predict,
                                  type = 'response'))

      weights <- 1/weights[resp_level]
            weights
    })

  .weighted_adjsfit(w, adjs, weight_lim, type, time)

}




ipw_rf <- function( adjs,formula = NULL, type ='surv', use_cont = TRUE, time = NULL, weight_lim = NULL, ...){

  w <- try({
      intrst <- Interest(adjs)

      if (use_cont){
        adjvars <- Adjvars(continous(adjs))

      } else {

        adjvars <- Adjvars(adjs)

      }

      if(is.null(formula)){
        nms <- names(adjvars)
        formula <- as.formula(paste0('1 ~ ', paste0(nms, collapse = ' + ')))
      }

      len <- length(intrst)
      resp_level <- lapply(levels(intrst),
                           function(x) x == intrst)
      resp_level <- lapply(resp_level, as.factor)
      resp_level <- unname(as.data.frame(resp_level))

      sampsize <- ceiling(rep(min(table(intrst)),2) * .632)

      def_params <- list(replace = FALSE, sampsize = sampsize)

      rf_params <- list(...)

      params <- modifyList(def_params, rf_params)

      model_level <- lapply(resp_level, function(x){

                              formula <- update.formula(formula, x ~ . )

                              lev_params <- list(formula = formula,
                                                 data = data.frame(x, adjvars))

                              params <- modifyList(lev_params, params)

                              do.call('randomForest', params)

                            }  )


      weights <- lapply(model_level, predict, type = 'prob')
      weights <- lapply(weights, function(x) x[ ,2])
      weights <- simplify2array(weights)

      resp_level <- lapply(resp_level, as.character)
      resp_level <- lapply(resp_level, as.logical)
      resp_level <- lapply(resp_level, as.numeric)
      resp_level <- simplify2array(resp_level)

      int_idx <- cbind(1:nrow(weights), apply(resp_level,1, which.max))
      weights <- 1/weights[int_idx]
      weights
    })

  .weighted_adjsfit(w, adjs, weight_lim, type, time)

}

## Resampling

.resample_if_different <- function(x, size){
  if(length(x) == size){
    x
  } else {
    x[sample(x = length(x), size = size, rep = TRUE) ]
  }
}


resampling <- function( adjs, ref = NULL, time = NULL, ...){

  s <- try({

      if(is.null(ref)){
        ref <- table(Adjvars(adjs))
      }

      sobj <- Sobj(adjs)
      covrts <- Covariates(adjs)
      intrst <- Interest(adjs)
      time <- time(adjs)

      covrts_groups <- split( 1:nrow(covrts), covrts)
      resampl_groups <- Map(.resample_if_different,
                            x = covrts_groups, size = ref)

      resampl <- unlist(resampl_groups, use.names = FALSE)
      resampl_sobj <- sobj[resampl, ]
      resampl_intrst <- intrst[resampl]

      sfit <- Sfit(resampl_sobj, resampl_intrst)
      s <- scurve_survfit(sfit, time)

    })

  .adjsfit(s)
}



## Conditional Methods

stratif <- function(adjs, time =  NULL, ...){

  s <- try({
  covrts <-  Covariates(adjs)
  adjvars <- Adjvars(adjs)
  sobj <- Sobj(adjs)
  if( is.null(time)) time <-  Sobj(adjs)[ ,'time' ]

  allstrat <- as.factor(as.integer(interaction(covrts)))

  alltab <- unique(cbind(inter = allstrat,covrts))
  adjtab <- data.frame(prop.table(table(adjvars)))

  tab <- merge(alltab, adjtab)

  sfit <- survfit(sobj ~ allstrat)
  inter <- as.integer(gsub('\\D',  '', names(sfit$strata)))
  mch <- match(x = inter, table = tab[ , c('inter')])

  sm <- summary(sfit, times = time, extend  = TRUE)

  tlen <- length(time)
  nstrat <- length(sfit$strata)
  p <- tab$Freq[mch]
  lev <- tab[mch,names(covrts)[1]]

  S <- sm$surv
  V <- sm$std.err^2
  dim(S) <- c(tlen, nstrat)
  dim(V) <- c(tlen, nstrat)

  w <- matrix(nrow = nstrat, ncol = nlevels(lev))
  w[ ,1:nlevels(lev)] <- p
  w <- w * model.matrix( ~ 0 + lev)

  s <- S %*% w
  se <- sqrt(V %*% w)

  group <- factor(rep(levels(lev), each = nrow(s)))

  s <- data.frame(time = rep(sort(time), ncol(s)),
                  surv = c(s),
                  se = c(se),
                  group = group)
  s
    })
  .adjsfit(s)
}

# Modelling

cox <- function(adjs, formula = NULL,  use_cont = TRUE, time = NULL, weight_lim = NULL, ...){

  s <- try({
      intrst <- Interest(adjs)
      sobj <- Sobj(adjs)
      if( is.null(time)) time <-  Sobj(adjs)[ ,'time' ]

      # using continous variable
      if (use_cont){
        covrts <- Covariates(continous(adjs))

      } else {

        covrts <- Covariates(adjs)

      }

      dt <-  cbind(sobj, covrts)
      #default formula
      if(is.null(formula)){
        nms <- names(Covariates(adjs))
        frm <- as.formula(paste0('1 ~ ', paste0(nms, collapse = ' + ')))
      } else {
        frm <- formula
      }

      strata <- ! is.null(attr(terms(frm, specials = 'strata'), 'specials')$strata) #detect stratified models




      frm <- update.formula(frm, sobj ~ .)


      coxfit <- coxph(formula(deparse(frm)), data = dt, singular.ok = TRUE )

      coxfit$call[2] <- call(deparse(frm))

      if(strata){

        s <- survfit(coxfit)
        s <- scurve_survfit(s, time)

      } else  {

        survxp <- survexp(  ~ intrst , ratetable = coxfit, data = dt,
                          times = sort(time))

        s <- scurve_survexp(survxp)
      }
      s
    })
  .adjsfit(s)

}




########################################
# Scenario definition


define_scen <- function(adjs, adjust, n = 10, ref = NULL, cont_vars =TRUE, weight_lim =c(.001, 1000), cens_type ='expected', max_t = NULL){

  if( n < 2) stop('In order to compare performance, a minimum of two replicas for scenario is required')

  adjvars <- Adjvars(continous(adjs))
  intrst <- Interest(adjs)

  sc <- paste0('sc_',levels(intrst))
  sh <- paste0('sh_',levels(intrst))
  mod <- model.matrix( ~ ., data = adjvars)[ ,-1 ,drop = FALSE]
  bet <- colnames(mod)


  nms <- c(sc,sh,bet,'cens','dataset')
  frame <- data.frame(t(NA[seq_along(nms)]), stringsAsFactors= TRUE)
  frame <- frame[-1, ]
  names(frame) <- nms

  init_list <- list(kaplan_meier = list() )
  adjust <- modifyList(init_list,adjmet)

  adjust <- lapply(adjust, modifyList,
                   list(adjs = substitute(adjs),
                        ref = ref,
                        weight_lim = weight_lim,
                        continous = cont_vars))

  ret <- list(scen = frame, adjust = adjust)

  attr(ret,'scl_col') <- seq_along(sc)
  attr(ret,'shp_col') <- length(sc) + seq_along(sh)
  attr(ret,'bet_col') <- length(c(sc,sh)) + seq_along(bet)
  attr(ret,'c_col') <- length(nms) - 1
  attr(ret,'dataset_col') <- length(nms)

  attr(ret,'adjsurv') <- substitute(adjs)
  attr(ret,'n') <- n
  attr(ret,'ref') <- ref
  attr(ret,'weight_lim') <- weight_lim
  attr(ret,'cont_vars') <- cont_vars
  attr(ret,'cens_type') <- cens_type
  attr(ret,'max_t') <- max_t
  attr(ret,'interest_lev') <- levels(intrst)

  class(ret) <- c(class(ret), 'simscen')

  ret
}

Adjs <- function(simscen){
  eval(attr(simscen,'adjsurv'))
}

time.simscen <- function(simscen){
  attr(simscen,'time')
}


print.simscen <- function(simscen){

  print(simscen$scen)

}

Scenarios <- function(simscen){
simscen$scen
}


cens_type <- function(simscen){
attr(simscen,'cens_type')
}

`cens_type<-` <- function(simscen, value){
  attr(simscen,'cens_type') <- value
  simscen
}



complete_scen <- function(simscen, limit = TRUE){

  scen <- Scenarios(simscen)
  lscen <- vapply(scen, length, integer(1))
  allfill <- all(lscen  > 0)

  if(allfill){

    if(limit){

      len <- prod(sapply(scen, length))

      if( len > 1000){

      stop('Combination of all parameters will produce ',len , ' scenarios which may be slow. Maybe you prefer shake(scenarios).\n If more than 1000 scenarios desired, turn off limit parameter')

    }

      simscen$scen <- expand.grid(scen,
                                  stringsAsFactors= FALSE)
      attr(simscen$scen,"out.attrs") <- NULL
      simscen
    }
  } else {

    stop('Can not be completed. There are missing parameters')

  }


}


betas <- function(scen){
  pos <- attr(scen,'bet_col')
  scen$scen[ ,pos]
}


`betas<-` <- function(x , value){
  pos <- attr(x,'bet_col')
  nparm <- length(pos)
  atr <- attributes(x)

  if ( length(value) == nparm ){

    x$scen <- as.list(x$scen)
    x$scen[pos] <- value
    len <- vapply(x$scen, function(x) length(x), integer(1))

    if( length(unique(len)) == 1 ){

      atr$row.names <- NULL
      scen  <- data.frame(x$scen ,stringsAsFactors= FALSE)
      x$scen <- scen
      cat('simulation parameters complete \n')
      class(x) <- c('simscen','list')
      x

    } else {

      class(x) <- c('list','simscen')
      x

    }

     } else {

       stop('list should be length ', nparm)

     }
}

scales <- function(scen){
  pos <- attr(scen,'scl_col')
  scen$scen[ ,pos]
}


`[.simscen` <- function(x, i, j){
  x$scen  <- x$scen[i , j]
  x
}




`scales<-` <- function(x , value){
   pos <- attr(x,'scl_col')
  nparm <- length(pos)
  atr <- attributes(x)


      if ( length(value) == nparm ){

    x$scen <- as.list(x$scen)
    x$scen[pos] <- value
    len <- vapply(x$scen, function(x) length(x), integer(1))

    if( length(unique(len)) == 1 ){

      atr$row.names <- NULL
      scen  <- data.frame(x$scen ,stringsAsFactors= FALSE)
      x$scen <- scen
      cat('simulation parameters complete \n')
      class(x) <- c('simscen','list')
      x

    } else {

      class(x) <- c('list','simscen')
      x
    }

     } else {

       stop('list should be length ', nparm)

     }
}



shapes <- function(scen){
  pos <- attr(scen,'shp_col')
  scen$scen[ ,pos]
}



`shapes<-` <- function(x , value){
  pos <- attr(x,'shp_col')
  nparm <- length(pos)
  atr <- attributes(x)

  if ( length(value) == nparm ){

    x$scen <- as.list(x$scen)
    x$scen[pos] <- value
    len <- vapply(x$scen, function(x) length(x), integer(1))

    if( length(unique(len)) == 1 ){

      atr$row.names <- NULL
      scen  <- data.frame(x$scen ,stringsAsFactors= FALSE)
      x$scen <- scen
      cat('simulation parameters complete \n')
      class(x) <- c('simscen','list')
      x

    } else {


      class(x) <- c('list','simscen')
      x

    }

     } else {

       stop('list should be length ', nparm)

     }
}


cens <- function(scen){
  pos <- attr(scen,'c_col')
  scen$scen[ ,pos]
}




`cens<-` <- function(x , value){
  pos <- attr(x,'c_col')
  nparm <- length(pos)
  atr <- attributes(x)

    if ( length(value) == nparm ){

    x$scen <- as.list(x$scen)
    x$scen[pos] <- value
    len <- vapply(x$scen, function(x) length(x), integer(1))

    if( length(unique(len)) == 1 ){

      atr$row.names <- NULL
      scen  <- data.frame(x$scen ,stringsAsFactors= FALSE)
      x$scen <- scen
      cat('simulation parameters complete \n')
      class(x) <- c('simscen','list')
      x

    } else {


      class(x) <- c('list','simscen')
      x

    }

     } else {

       stop('list should be length ', nparm)

     }
}

dataset <- function(scen){
  pos <- attr(scen,'dataset_col')
  scen$scen[ ,pos]
}





`dataset<-` <- function(x , value){
  pos <- attr(x,'dataset_col')
  nparm <- length(pos)
  atr <- attributes(x)

   if ( length(value) == nparm ){

    x$scen <- as.list(x$scen)
    x$scen[pos] <- value
    len <- vapply(x$scen, function(x) length(x), integer(1))

    if( length(unique(len)) == 1 ){

      atr$row.names <- NULL
      x$scen  <- data.frame(x$scen ,stringsAsFactors= FALSE)

      cat('simulation parameters complete \n')
      class(x) <- c('simscen','list')
      x


    } else {


      class(x) <- c('list','simscen')
      x

    }

     } else {

       stop('list should be length ', nparm)

     }
}


#Simulate function
########################################



weib_sim <- function(simscen){

  scl <- t(scales(simscen))
  shp <- t(shapes(simscen))
  bet <- t(betas(simscen))
  cns <- as.vector(cens(simscen))
  cens_type <- attr(simscen,'cens_type')
  max_t <- attr(simscen,'max_t')

  data <- get(dataset(simscen))

  adjs <- AdjSurv(formula(Adjs(simscen)), NULL , data)
  adjvars <- Adjvars(continous(adjs))
  intrst <-  Interest(adjs)

  mf <- model.frame(~. ,adjvars)
  X <- model.matrix(mf, adjvars)[ ,-1, drop = FALSE]

  weib_param <- data.frame(scl,shp)[intrst, ]
  names(weib_param) <- c('scale','shape')

  esc_lmbd <- weib_param$scale *
    exp(-( X %*% bet) /weib_param$shape)

  Tm <- Map(rweibull, scale = esc_lmbd,
           shape = weib_param$shape, n = 1)

  Tm <- unlist(Tm)

  l <- quantile(Tm, .632)
  k <- log(log(2))/(log(median(Tm))-log(l))

  if(is.null(cens_type)){

   # No censoring
    obs_T <- Tm
    delta <- rep(1, length(Tm))

  } else {
     # Random
    if(cens_type == 'expected'){

      opt <- optimize(
          function(x) {
            abs(mean(pweibull(Tm, scale = x, shape = k) - cns * 1.01))
          }, interval = c(1e-16,1e+16))
      l_cens <- opt$minimum
      k_cens <- k

    } else if (cens_type == 'factor'){
      l_cens <- l * cns[1]

      if( is.na( cns[2]) ) cns[2] <- 1
      k_cens <- k * cns[2]

    } else if (cens_type == 'random') {

      l_cens <- l * exp(rnorm(1))
      k_cens <- k * exp(rnorm(1))

    } else if (cens_type == 'informative') {

      l_cens <- esc_lmbd *  ((1-cns)^ (1 /weib_param$shape ))/
        (cns^ (1 /weib_param$shape ))
      k_cens <- weib_param$shape

    }

    cens_T <- rweibull(length(Tm), scale = l_cens , shape = k_cens)
    delta <- Tm < cens_T
    obs_T <- ifelse(delta, Tm, cens_T)
  }

  not_observed <- obs_T > max_t
  obs_T[not_observed] <- max_t
  delta[not_observed] <- 0

  delta <- as.numeric(delta)

  adjs$Surv <- Surv(obs_T, delta)
  adjs
}


# Reference functions
########################################
dweibull_scmix <- function(x, w, scales, shape, log = FALSE){
  if( length(w) == length(scales)){
  len <- length(x)
  res <- vapply(scales ,dweibull, numeric(len), x = x, shape = shape, log = log) %*% w
  c(res)
} else {
  stop(' w not equal scale parameters')
}
}



pweibull_scmix <- function(x, w, scales, shape, log = FALSE){
  if( length(w) == length(scales)){
  len <- length(x)
  res <- vapply(scales ,pweibull, numeric(len), q = x, shape = shape, log.p = log) %*% w
  c(res)
  } else {
  stop(' w not equal scale parameters')
}
}


sweibull_scmix <- function(x, w, scales, shape){
  if( length(w) == length(scales)){
  len <- length(x)
  res <- 1 - vapply(scales ,pweibull, numeric(len), q = x, shape = shape) %*% c(w)
  c(res)
    } else {
  stop(' w not equal scale parameters')
}
}


sweibfun_scmix <- function(w, scales, shape){
  S <- function(x){
    sweibull_scmix(x, w, scales = scales, shape = shape)
  }
    attr(S,'w') <- w
    attr(S,'scales') <- scales
    attr(S,'shape') <- shape
  S
}

weib_scurve <- function(simscen){

  scl <- scales(simscen)
  shp <- shapes(simscen)
  bet <- t(betas(simscen))

  data <- get(as.character(dataset(simscen)))

  adjs <- AdjSurv(formula(Adjs(simscen)),c(), data)

  adjvars <- Adjvars(continous(adjs))
  intrst <- Interest(adjs)
  nlev <- vapply(adjvars, nlevels, numeric(1))

  adjvars_cont <- adjvars[nlev < 1]

  ref_val <- vapply(adjvars_cont, function(x) {
                      mean(sample(x, 1500, replace = TRUE))
                    }
                  , double(1))

  adjvars_fact <- adjvars[nlev > 0]

  idx <- lapply(nlev[nlev > 0], function(x) 1:x)
  idx_cmb <- expand.grid(idx)
  X <- Map( function(x, nl)  diag(nl)[x, -1, drop =FALSE],
           x = idx_cmb, nl = nlev[ nlev > 0 ])
  X <- do.call('cbind',X)
  X <- cbind(ref_val, X)

  scaled_l <- Map(function(sc, sh) c(sc *
                                       exp(X %*% bet) ^(-1 /sh)),
                  sc = scl,
                  sh = shp)

  p_intrst <- list(prop.table(table(unique(adjvars_fact))))

  ret <- Map(sweibfun_scmix, p_intrst, scale =scaled_l, shape = shp)

  names(ret) <- levels(intrst)

  ret
}



Adjustings <- function(simscen){
  simscen$adjust
}

N <- function(simscen){
   attr(simscen, 'n')
}





.adjsim_row <- function(simscen, n, ret_sobj = FALSE){
  # adjs <- Adjs(simscen) # do not remove
  adjust <- Adjustings(simscen)
  nadjust <- length(adjust)

  scenpb <- txtProgressBar(0,1, style = 2)
  res <- lapply( seq_len(n), function(i){

                     sadjs <- weib_sim(simscen)
                     stime <- sort(time(sadjs))

                     calls <- lapply(adjust, modifyList,
                                     val = list(adjs = substitute(sadjs)))
                     sadjfit <- Map(do.call, what = names(calls), calls)


                     sfit <- survfit(Sobj(sadjs) ~  Interest(sadjs))
                     failed <- vapply(sadjfit, inherits, logical(1), 'try-error')

                     surv <- lapply(sadjfit, function(x) x$surv )
                     nrw <- length(stime) *  nlevels(Interest(sadjs))
                     curves <- matrix(NA, nrw, nadjust)
                     curves <- data.frame(curves)

                     wrong <- vapply(surv, length, numeric(1)) < nrw
                     failed[wrong] <- TRUE

                     curves[!failed] <- surv[!failed]
                     curves <- simplify2array(curves)

                     forms <- lapply(adjust, function(x)
                       if(is.null(x$formula)){
                         ''
                       } else {
                         paste0(':',as.character(x$formula[2]))

                       }    )

                     colnames(curves) <- Map(paste0,
                                            names(adjust),
                                            forms)


                     ref_fun <- weib_scurve(simscen)

                     ref_curves <- lapply(ref_fun, function(ref_f){
                                            ref_f(stime)
                                          }    )
                     ref_curves <- unlist(ref_curves)
                     ref_curves <- unname(ref_curves)
                     ref_curves <- matrix(ref_curves)[ ,rep(1,length(adjust))]

                     if(ret_sobj){
                       ret <- list(scurves = curves,
                                   reference =ref_curves,
                                   nrisk = summary(sfit,
                                       times = stime,
                                       extend = TRUE)$n.risk,
                                   time = summary(sfit,
                                       times = stime,
                                       extend = TRUE)$time,
                                   sobj = cbind(Sobj(sadjs),
                                       group = Interest(sadjs)),
                                   failed =failed)
                     } else {
                       ret <- list(scurves = curves,
                                   reference =ref_curves,
                                   nrisk = summary(sfit,
                                       times = stime,
                                       extend = TRUE)$n.risk,
                                   time = summary(sfit,
                                       times = stime,
                                       extend = TRUE)$time,
                                   failed =failed)
                     }
                     setTxtProgressBar(scenpb, i/ n)
                     ret

                   })

  res <- do.call('cbind',res)
  res <- apply(res,1,  simplify2array)

  res$failed <- t(res$failed)
  close(scenpb)
  res
}



shake <- function(simscen, n = 1 , what ='all' ,amount = 1/10){
  atr <- attributes(simscen)

  non_negative  <- c(atr$scl_col,
                     atr$shp_col,
                     atr$c_col)

  scen <- Scenarios(simscen)

  sel <- seq_along(scen)
  sel <- switch(what,
         all = sel[with(atr, c(shp_col, scl_col, bet_col, c_col))],
         shapes = sel[ atr$shp_col],
         scales = sel[ atr$scl_col],
         betas = sel[ atr$bet_col])


  shk_param <- replicate({

      psd <- amount * colMeans(abs(scen[ ,sel]))
      shk <- t(scen[sel]) + rnorm(n = ncol(scen[sel]) ,
                                               mean = 0, sd = psd)
      shk <- t(shk)
      idx <- sel %in% non_negative

      shk[ ,idx] <- abs(shk[ ,idx])

      scen[sel] <- shk

      scen$cens[scen$cens >= 1] <-
        c( 1 / scen$cens[scen$cens >= 1])

      scen

    }, n = n, simplify = FALSE)

  shk_param <- do.call('rbind',shk_param)

  simscen$scen <- rbind(simscen$scen, shk_param)
  simscen
}



AdjSim <- function(simscen, ret_sobj = FALSE, mc.cores = mc.control ){

  nscen <- nrow(Scenarios(simscen))
  nsim <- N(simscen)
  pb <- txtProgressBar(0,1, style = 3)
  cat('\n')
  adjsim <- mclapply(seq_len(nscen),
                   function(i){
                     ret <- .adjsim_row(simscen[i, ],
                                      n = nsim,
                                        ret_sobj = ret_sobj)
                     setTxtProgressBar(pb, i/nscen)
                     ret
                   }, mc.cores = mc.cores)

  attr(adjsim, "adjsurv") <-   attr(simscen, 'adjsurv')
  attr(adjsim, 'n') <-   attr( simscen, 'n')
  attr(adjsim,'interest_lev') <-   attr(simscen, 'interest_lev')
  attr(  adjsim, 'weight_lim') <-   attr(simscen, 'weight_lim')
  attr(  adjsim, 'cont_vars') <-   attr( simscen, 'cont_vars')
  attr( adjsim, 'cens_type') <-   attr(simscen, 'cens_type')
  attr(adjsim,'max_t') <-   attr(simscen, 'max_t')
  attr(adjsim, 'simscen') <-   simscen
  class(adjsim) <- c('adjsim')
  setTxtProgressBar(pb, 1)
  close(pb)

  attr(adjsim,'print') <- summary(adjsim)
  adjsim
}






Reference <- function(adjsim){
  lapply(adjsim, function(x) x$reference)
}


Curves <- function(adjsim){
  lapply(adjsim, function(x) x$scurves)
}

Sobj.adjsim <- function(adjsim){
  lapply(adjsim, function(x) x$sobj)
}


Failed <- function(adjsim, scen = TRUE){

  res <- lapply(adjsim, function(x) x$failed)

  if(scen){
   return(which( vapply(res, any, logical(1))))
  }
  simplify2array(res)
}


Nrisk <- function(adjsim){
  lapply(adjsim, function(x) x$nrisk)
}

levels.adjsim <- function(adjsim){
  attr(adjsim, 'interest_lev')
}

time.adjsim <- function(adjsim){
  lapply(adjsim, function(x) x$time)
}





`[.adjsim` <- function(x, i, j, k, mc.cores = mc.control){

  mi <- missing(i)
  mj <- missing(j)
  mk <- missing(k)

  atr <- attributes(x)

  if(mi) i <- seq_along(x)
  ret <- lapply(i, function(idx) x[[idx]])

  if(mj) j <- seq_len( ncol(ret[[1]][[1]]))
  if(mk) k <- seq_len(dim(ret[[1]][[1]])[3])

  ret <- mclapply(ret, function(obj){
                  list(scurves = obj$scurves[ ,j,k, drop = FALSE],
                       reference = obj$reference[ ,j, k, drop = FALSE],
                       nrisk = obj$nrisk[ ,  k, drop = FALSE],
                       time = obj$time[ ,k, drop = FALSE],
                       sobj = obj$sobj,
                       failed = obj$failed[k, j, drop = FALSE])
  }, mc.cores = mc.cores    )

  attributes(ret) <- atr
  attr(ret, 'print') <- NULL
  ret
}



residuals.adjsim <- function(adjsim, squared =TRUE, nrisk_w = NULL, mc.cores = mc.control){
  sq <- 2 ^squared

  ret <- mcMap( function(x,y) (x - y )^sq ,
               x  = Reference(adjsim),
               y = Curves(adjsim), mc.cores = mc.cores)

  if( is.null(nrisk_w)){

    return(ret)

  } else {

    if(inherits(nrisk_w, "function")){

      wn <- lapply(Nrisk(adjsim), nrisk_w)

    } else {

      wn <- Nrisk(adjsim)
    }


   ret <- mcMap(function(ajs, wl){
          lapply( seq_len(ncol(wl)), function(i){
                   ajs[ , ,i] * wl[ ,i]
                 }    )
        },
              ajs = ret,
              wl = wn, mc.cores = mc.cores)

    lapply(ret, simplify2array)

  }
}



print.adjsim <- function(adjsim){

   cat('\n')

    cat('\n MSE after ')
    cat(N(adjsim))
    cat(' simulations of ')
    cat(length(adjsim))
   cat(' scenarios: \n \n ')

  if(is.null(attr(adjsim, 'print'))){
    print(summary(adjsim))

  } else {

    print(attr(adjsim,'print'))
    cat('\n')
  }

}



plot.adjsim <- function(adjsim, which , scen_summ = mean, global_summ = mean, nrisk_w = NULL, squared = TRUE, plot = TRUE, mc.cores = mc.control, ... ){

  scen <- mclapply( resid(adjsim, squared= squared, nrisk_w = nrisk_w),
                 function(ajs){
                   apply(ajs,1:2  , scen_summ)
                 }, mc.cores = mc.cores )

  scen <- simplify2array(scen)


  scen <- apply(scen,1:2, global_summ)

  group <- factor(rep(levels(adjsim), each = nrow(scen) / nlevels(adjsim)))

  if(missing(which)) which <- seq_along(levels(adjsim))
  wlen <- length(which)

  if(wlen > 2 && wlen %% 2 > 0) {
    lo <- c(seq_along(which),0)
    lo <- matrix( lo , ncol = 2, byrow = TRUE)
  } else {
    lo <- seq_along(which)
    lo <- matrix(lo, ncol = 2, byrow = TRUE)
  }

  ol <- layout(lo)
  op <- par( cex =.7)
  invisible({
  lapply(levels(adjsim)[which], function(lvl){

           if(squared) yl <- 'Averaged square error'
           else yl <- 'Averaged error'

           xl <- 'Ordered observed times'

           matplot(scen[group == lvl, ], type ='l',
                   ylab = yl,
                   xlab = xl,
                   col = 1:ncol(scen),
                   lty = 1:ncol(scen),
                   ...)
           abline(h =0, col ='grey', lty =3)
           title(lvl)
           legend('topright', colnames(scen),lty =1:ncol(scen),
                  col = 1:ncol(scen), cex =.9, bty ='n')
         }    )
    })
  par(op)
  layout(1)
}


survplot <- function(adjsim, scenario = 1, replica = 1, adjust, which, ...){

  if(missing(adjust)) adjust <- seq_len(dim(adjsim[[1]][[1]])[2])
  subset <- adjsim[scenario, adjust, replica]
  scurves <- Curves(subset)[[1]]
  nms <- colnames(scurves)
  dim(scurves) <- dim(scurves)[1:2]

  ref <- Reference(subset)[[1]][ ,1, ]
  dim(ref) <- dim(ref)[1:2]
  scurves <- cbind(ref,scurves)
  colnames(scurves) <- c('reference', nms)

  group <- factor(rep(levels(adjsim),
                      each = nrow(scurves) / nlevels(adjsim)))


  if(missing(which)) which <- seq_along(levels(adjsim))
  wlen <- length(which)

  if(wlen > 2 && wlen %% 2 > 0) {
    lo <- c(seq_along(which),0)
    lo <- matrix( lo , ncol = 2, byrow = TRUE)
  } else {
    lo <- seq_along(which)
    lo <- matrix(lo, ncol = 2, byrow = TRUE)
  }


  ol <- layout(lo)
  op <- par(cex =.7)
  invisible({
      sq_lv <- seq_len(nlevels(adjsim)+1)
  lapply(levels(adjsim)[which], function(lvl){
           xl <- 'Ordered observed Time'
           yl <- 'Survival'
           matplot(scurves[group == lvl, ], type ='s',
                   ylab = yl,
                   xlab = xl, ...)
           title(lvl)
           legend('topright',
                  colnames(scurves),
                  lty =sq_lv,
                  col = sq_lv,
                  cex =.8,
                  bty ='n')
         }    )
    })
  par(op)
  layout(1)

}



#Summary functions

test_scen <- function(sm){

  sm <- t(sm)

  fritest <- friedman.test(sm)$p.value

  ncl <- ncol(sm)
  sq_met <- seq(ncl)[-1]

  nas <- apply(apply(sm, 2, is.na), 2 , all)[-1]
  fill <- sq_met[!nas] -1

  nas <- which(nas) + 1

  ttest <- NA[sq_met]
  wtest <- ttest

  ttest[fill] <- apply(sm[  ,- c(1, nas), drop = FALSE] - sm[ ,1], 2, function(x) t.test(x)$p.value)

  names(ttest) <- paste0('tt.',colnames(sm)[-1])

  wtest[fill] <- apply(sm[ ,-c(1,nas), drop = FALSE] - sm[ ,1], 2, function(x) wilcox.test(x)$p.value)

  names(wtest) <- paste0('wcx.',colnames(sm)[-1])

  mrank <- rowMeans(apply(sm,1,rank))

  names(mrank) <- paste0('rank.',names(mrank))

  list(friedman.test = fritest, ttest= ttest, wilcoxon= wtest,mranks = mrank)

}


## summ_test <- function(x, signif){
##   nr <- nrow(x)
##   nmet <- (nr - 2) / 3
##   ntest <- 1 + nmet * 2

##   sign <- apply(x[1:ntest, ,drop = FALSE ] <  signif,1,sum)/
##     ncol(x)

##   rnk <- rowMeans(x[(ntest+1) :nr, , drop = FALSE])
##   c(sign, rnk)
## }

## summ_test <- function(x, signif){



##   rnk <- rowMeans(x[(ntest+1) :nr, , drop = FALSE])
##   c(sign, rnk)
## }


summary.adjsim <- function(adjsim, rm.failed = FALSE, squared = TRUE, nrisk_w = NULL, conf = 0.95, mc.cores = mc.control){

  ret <- .adjsimsm_scen(adjsim,
                        squared = squared,
                        nrisk_w =nrisk_w,
                        rm.failed = rm.failed,
                        mc.cores = mc.cores)


  failed <- Failed(adjsim , scen = FALSE)
  wf <- unique(which(failed, arr.ind = TRUE)[ ,2])

  if( length(wf) > 0 ){
    fl <- colnames(failed)[wf]

    text <- paste0('\n\nFollowing methods have failed scenarios:\n',
                   paste(colnames(failed)[wf], collapse= ', '), '\n',
                   'Summary is computed with available simulations.')
    warning(text)

  }
  ret
}



.adjsimsm <- function(adjsim, squared, nrisk_w, conf, mc.cores){
  try({
  cat('summarizing:\n')
  pb <- txtProgressBar(0,1, style =3)
  resd <- resid(adjsim, nrisk_w = nrisk_w,squared = squared)

  setTxtProgressBar(pb, 1/5)
  mean_scen <- mclapply(resd, function(x){
                          nms <- colnames(x)
                          d <- dim(x)
                          ret <- colMeans(x, na.rm = TRUE,
                                          dims= 1)
                          rownames(ret) <- nms
                          ret
                        }            , mc.cores = mc.cores)

  setTxtProgressBar(pb, 2/5)

  sm <- simplify2array(mean_scen)
  nms <- rownames(sm)
  d <- dim(sm)
  nwd <- c(d[1], d[2] * d[3])
  dim(sm) <- nwd
  me <- rowMeans(sm, na.rm = TRUE)

  setTxtProgressBar(pb, 3/5)

  se <- apply(sm, 1, sd)
  sm <- rbind(me, se)
  colnames(sm) <- nms
  rownames(sm) <- c('mean', 'se')

  setTxtProgressBar(pb, 4/5)

  test <- mclapply(mean_scen, test_scen, mc.cores = mc.cores)
  test <- do.call('cbind',test)
  test <- apply(test,1, simplify2array)
  ranksm <- rowMeans(test[[4]])

  test[[1]] <- matrix(test[[1]], nrow = 1)
  sign <- lapply(test[1:3], function(x){
                   rowMeans(x  < (1 - conf), na.rm = TRUE)
                 } )

  setTxtProgressBar(pb, 5/5)
  close(pb)
      c(sm =list(sm), sign, ranksm =list(ranksm))
    })
}



.adjsimsm_scen <- function(adjsim, groups = NULL , squared, nrisk_w, rm.failed, mc.cores){

  try({
      cat('\n summarizing:\n')
      pb <- txtProgressBar(0,1, style =3)

      if(rm.failed){
        failed <- Failed(adjsim)
        keep <- setdiff(seq(length(adjsim)), failed)
        adjsim <- adjsim[keep]
      }

      setTxtProgressBar(pb, 1/5)
      resd <- resid(adjsim, nrisk_w = nrisk_w,squared = squared)


      setTxtProgressBar(pb, 2/5)

      if(is.null(groups)) groups <- seq_len(length(adjsim))
      scen_mean <- mclapply(resd, colMeans , mc.cores = mc.cores,
                            na.rm = TRUE, dim = 1)

      setTxtProgressBar(pb, 3/5)
      sm <- mclapply(scen_mean, rowMeans, mc.cores = mc.cores)
      sm <- simplify2array(sm)

      setTxtProgressBar(pb, 4/5)
      test <- mclapply(scen_mean, test_scen, mc.cores = mc.cores)
      test <- simplify2array(test)
      test <- apply(test, 1, simplify2array)

      setTxtProgressBar(pb, 5/5)
      close(pb)
      cat('\n')
      ret <- list(sm = sm, comp = test)
      class(ret) <- c('adjssm',class(ret))
      ret
    })
}

print.adjssm <- function(x, maxprint = 5){

  nc <- ncol(x$sm)
  nprint <- seq_len(min(nc, maxprint))

  print(  list(error = x$sm[ ,nprint],
               friedman.test =x$comp$friedman.test[nprint],
               t.test = x$comp$ttest[ ,nprint],
               wilcoxon= x$comp$wilcoxon[ ,nprint],
               mean.ranks =x$comp$mranks[ ,nprint]))

  if(nc > maxprint){
    txt <- paste('...Printing',
                 maxprint,
                 'out of', nc, 'scenario groups \n')
    cat(txt)
      }
}


plot.adjssm <- function(adjssm, g = NULL, ranked = FALSE, conf = .95, ...){

  if(ranked){

    rnk <- apply(adjssm$comp$mranks, 2, rank)

  } else {

    rnk <- adjssm$comp$mranks

  }

  nr <- nrow(rnk)
  nc <- ncol(rnk)

  if( is.null(g) ) g <- seq_len(ncol(rnk))

  friedm <- rep(NA,nc)
  friedm[which(adjssm$comp$friedman.test < (1-conf))] <- 8


  wlcxs <- which(adjssm$comp$wilcoxon < (1 - conf))
  tts <- which(adjssm$comp$ttest < (1 - conf))

  pchidx <- rnk[-1, ]
  pchidx[] <- NA

  pchidx[wlcxs] <- 1
  pchidx[tts] <- 3
  pchidx[wlcxs[wlcxs %in% tts]] <- 10

  reset_par <- par(mar = c(5, 4, 4, 6))
  matplot(y = t(rnk), x = g, type = 'l',
          ylim =c(nr ,1), bty = 'L')

invisible({
  Map(function(adjssm, pch, cl){
        points(y = adjssm, x =g, pch = pch, col = cl)
      } ,
      y = data.frame(t(rnk[-1, ])),
      pch = data.frame(t(pchidx)),
      cl = c(2:6,1 : (nr -6)))
  })

  points(x = g, y = rep( nr - .1, length(g)),
         pch = friedm,
         cex =.8)

  legend(x = nc+ nc * 1e-2 /2, y = 1, rownames(y$sm),
         title = 'Adjusting method',
         xpd = TRUE,
         cex = .7,
         col = 1:6,
         lty = 1:5,
         box.lwd = 0)

  legend(x =nc + nc * 1e-2, y = nr/2 -.5 ,c('wilcoxon', 'T test', 'Friedman'),
         xpd = TRUE,
         cex = .7,
         pch = c(1,3, 8),
         box.lwd = 0, title = paste('p-val <', (1 - conf)))

  par(reset_par)


}

