mle3 <- function (minuslogl, start, method, optimizer, fixed = NULL,
          data = NULL, subset = NULL, default.start = TRUE, eval.only = FALSE,
          vecpar = FALSE, parameters = NULL, parnames = NULL, skip.hessian = FALSE,
          hessian.opts = NULL, use.ginv = TRUE, trace = FALSE, browse_obj = FALSE,
          gr = NULL, optimfun, ...)
{
  if (missing(method))
    method <- mle2.options("optim.method")
  if (missing(optimizer))
    optimizer <- mle2.options("optimizer")
  L <- list(...)
  if (optimizer == "optimize" && (is.null(L$lower) || is.null(L$upper)))
    stop("lower and upper bounds must be specified when using\n'optimize'")
  if (inherits(minuslogl, "formula")) {
    pf <- function(f) {
      if (is.null(f)) {
        ""
      }
      else {
        paste(f[2], "~", gsub(" ", "", as.character(f[3])),
              sep = "")
      }
    }
    if (missing(parameters)) {
      formula <- pf(minuslogl)
    }
    else {
      formula <- paste(pf(minuslogl), paste(sapply(parameters,
                                                   pf), collapse = ", "), sep = ": ")
    }
    tmp <- calc_mle2_function(minuslogl, parameters, start = start,
                              parnames = parnames, data = data, trace = trace)
    minuslogl <- tmp$fn
    start <- tmp$start
    fdata <- tmp$fdata
    parameters <- tmp$parameters
  }
  else {
    formula <- ""
    fdata <- NULL
  }
  call <- match.call()
  call.orig <- call
  call$data <- eval.parent(call$data)
  call$upper <- eval.parent(call$upper)
  call$lower <- eval.parent(call$lower)
  call$gr <- eval.parent(call$gr)
  call$control <- eval.parent(call$control)
  call$method <- eval.parent(call$method)
  if (!missing(start))
    if (!is.list(start)) {
      if (is.null(names(start)) || !is.vector(start))
        stop("'start' must be a named vector or named list")
      vecpar <- call$vecpar <- TRUE
      start <- as.list(start)
    }
  if (missing(start) && default.start)
    start <- formals(minuslogl)
  if (!is.null(fixed) && !is.list(fixed)) {
    if (is.null(names(fixed)) || !is.vector(fixed))
      stop("'fixed' must be a named vector or named list")
    fixed <- as.list(fixed)
  }
  if (!is.null(data) && !is.list(data))
    stop("'data' must be a list")
  nfix <- names(unlist(bbmle::namedrop(fixed)))
  if (!is.null(parnames(minuslogl))) {
    nfull <- parnames(minuslogl)
    fullcoef <- vector("list", length(nfull))
    names(fullcoef) <- nfull
  }
  else {
    fullcoef <- formals(minuslogl)
    nfull <- names(fullcoef)
  }
  if (any(!nfix %in% nfull))
    stop("some named arguments in 'fixed' are not arguments to the specified log-likelihood function")
  if (length(nfix) > 0)
    start[nfix] <- NULL
  fullcoef[nfix] <- fixed
  nstart <- names(unlist(sapply(bbmle::namedrop(start), eval.parent)))
  fullcoef[!nfull %in% nfix & !nfull %in% nstart] <- NULL
  nfull <- names(fullcoef)
  lc <- length(call$lower)
  lu <- length(call$upper)
  npnfix <- sum(!nfull %in% nfix)
  if (!npnfix == 0 && (lu > npnfix || lc > npnfix)) {
    warning("length mismatch between lower/upper ", "and number of non-fixed parameters: ",
            "# lower=", lc, ", # upper=", lu, ", # non-fixed=",
            npnfix)
  }
  template <- lapply(start, eval.parent)
  if (vecpar)
    template <- unlist(template)
  start <- sapply(bbmle::namedrop(start), eval.parent)
  nstart <- names(unlist(bbmle::namedrop(start)))
  oo <- match(nstart, names(fullcoef))
  if (any(is.na(oo)))
    stop("some named arguments in 'start' are not arguments to the specified log-likelihood function")
  start <- start[order(oo)]
  fix_order <- function(c1, name, default = NULL) {
    if (!is.null(c1)) {
      if (length(unique(c1)) > 1) {
        if (is.null(names(c1)) && length(unique(c1)) >
            1) {
          warning(name, " not named: rearranging to match 'start'")
          oo2 <- oo
        }
        else oo2 <- match(names(unlist(bbmle::namedrop(c1))),
                          names(fullcoef))
        c1 <- c1[order(oo2)]
      }
    }
    else c1 <- default
    c1
  }
  call$lower <- fix_order(call$lower, "lower bounds", -Inf)
  call$upper <- fix_order(call$upper, "upper bounds", Inf)
  call$control$parscale <- fix_order(call$control$parscale,
                                     "parscale")
  call$control$ndeps <- fix_order(call$control$ndeps, "ndeps")
  if (is.null(call$control))
    call$control <- list()
  denv <- local(environment(), c(as.list(data), fdata, list(mleenvset = TRUE)))
  argnames.in.data <- names(data)[names(data) %in% names(formals(minuslogl))]
  args.in.data <- lapply(argnames.in.data, get, env = denv)
  names(args.in.data) <- argnames.in.data
  args.in.data
  objectivefunction <- function(p) {
    if (browse_obj)
      browser()
    l <- bbmle::relist2(p, template)
    names(p) <- nstart[order(oo)]
    l[nfix] <- fixed
    if (vecpar) {
      l <- bbmle::namedrop(l[nfull])
      l <- unlist(l)
      args <- list(l)
      args <- c(list(l), args.in.data)
    }
    else {
      args <- c(l, args.in.data)
    }
    do.call("minuslogl", bbmle::namedrop(args))
  }
  objectivefunctiongr <- if (!is.null(gr))
    function(p) {
      if (browse_obj)
        browser()
      l <- bbmle::relist2(p, template)
      names(p) <- nstart[order(oo)]
      l[nfix] <- fixed
      if (vecpar) {
        l <- bbmle::namedrop(l[nfull])
        l <- unlist(l)
        args <- list(l)
        args <- c(list(l), args.in.data)
      }
      else {
        args <- c(l, args.in.data)
      }
      v <- do.call("gr", args)
      if (is.null(names(v))) {
        if (length(v) == length(l) && !is.null(tt <- names(l))) {
          vnames <- tt
        }
        else if (length(v) == length(p) && !is.null(tt <- names(p))) {
          vnames <- tt
        }
        else if (!is.null(tt <- parnames(minuslogl))) {
          vnames <- tt
        }
        else vnames <- names(formals(minuslogl))
        if (length(vnames) != length(v))
          stop("name/length mismatch in gradient function")
        names(v) <- vnames
      }
      v[!names(v) %in% nfix]
    }
  if (!("mleenvset" %in% ls(envir = environment(minuslogl)))) {
    newenv <- new.env(hash = TRUE, parent = environment(minuslogl))
    d <- as.list(denv)
    mapply(assign, names(d), d, MoreArgs = list(envir = newenv))
    environment(minuslogl) <- newenv
    if (!is.null(gr)) {
      newenvgr <- new.env(hash = TRUE, parent = environment(minuslogl))
      mapply(assign, names(d), d, MoreArgs = list(envir = newenvgr))
      environment(gr) <- newenvgr
    }
  }
  if (length(start) == 0 || eval.only) {
    if (length(start) == 0)
      start <- numeric(0)
    optimizer <- "none"
    skip.hessian <- TRUE
    oout <- list(par = start, value = objectivefunction(start),
                 hessian = matrix(NA, nrow = length(start), ncol = length(start)),
                 convergence = 0)
  }
  else {
    oout <- switch(optimizer, optim = {
      arglist <- list(...)
      arglist$lower <- arglist$upper <- arglist$control <- NULL
      do.call("optim", c(list(par = start, fn = objectivefunction,
                              method = method, hessian = FALSE, gr = objectivefunctiongr,
                              control = call$control, lower = call$lower, upper = call$upper),
                         arglist))
    }, optimx = {
      arglist <- list(...)
      arglist$lower <- arglist$upper <- arglist$control <- NULL
      do.call("optimx", c(list(par = start, fn = objectivefunction,
                               method = method, hessian = FALSE, gr = objectivefunctiongr,
                               control = call$control, lower = call$lower, upper = call$upper),
                          arglist))
    }, nlm = nlm(f = objectivefunction, p = start, hessian = FALSE,
                 ...), nlminb = nlminb(start = start, objective = objectivefunction,
                                       hessian = NULL, ...), constrOptim = constrOptim(theta = start,
                                                                                       f = objectivefunction, method = method, ...), optimize = ,
    optimise = optimize(f = objectivefunction, interval = c(call$lower,
                                                            call$upper), ...), user = {
                                                              arglist <- list(...)
                                                              arglist$lower <- arglist$upper <- arglist$control <- NULL
                                                              do.call(optimfun, c(list(par = start, fn = objectivefunction,
                                                                                       method = method, hessian = FALSE, gr = objectivefunctiongr,
                                                                                       control = call$control, lower = call$lower,
                                                                                       upper = call$upper), arglist))
                                                            }, stop("unknown optimizer (choices are 'optim', 'nlm', 'nlminb', 'constrOptim', 'user', and 'optimi[sz]e')"))
  }
  optimval <- switch(optimizer, optim = , constrOptim = , optimx = ,
                     user = , none = "value", nlm = "minimum", optimize = ,
                     optimise = , nlminb = "objective")
  if (optimizer == "optimx") {
    fvals <- oout[["value"]]
    conv <- oout[["convcode"]]
    best <- which.min(fvals)
    oout <- list(par = as.numeric(unlist(oout[best, 1:attr(oout,
                                                           "npar")])), value = fvals[best], convergence = conv[best],
                 method.used = attr(oout, "details")[, "method"][[best]])
  }
  if (optimizer == "nlm") {
    oout$par <- oout$estimate
    oout$convergence <- oout$code
  }
  if (optimizer %in% c("optimise", "optimize")) {
    oout$par <- oout$minimum
    oout$convergence <- 0
  }
  if (optimizer %in% c("nlminb", "optimise", "optimize") ||
      is.null(names(oout$par))) {
    names(oout$par) <- names(start)
  }
  if (length(oout$par) == 0)
    skip.hessian <- TRUE
  if (!skip.hessian) {
    if ((!is.null(call$upper) || !is.null(call$lower)) &&
        any(oout$par == call$upper) || any(oout$par == call$lower))
      warning("some parameters are on the boundary: variance-covariance calculations based on Hessian may be unreliable")
  }
  namatrix <- matrix(NA, nrow = length(start), ncol = length(start))
  if (!skip.hessian) {
    psc <- call$control$parscale
    if (is.null(gr)) {
      if (is.null(psc)) {
        oout$hessian <- try(hessian(objectivefunction,
                                    oout$par, method.args = hessian.opts))
      }
      else {
        tmpf <- function(x) {
          objectivefunction(x * psc)
        }
        oout$hessian <- try(hessian(tmpf, oout$par/psc,
                                    method.args = hessian.opts))/outer(psc, psc)
      }
    }
    else {
      if (is.null(psc)) {
        oout$hessian <- try(jacobian(objectivefunctiongr,
                                     oout$par, method.args = hessian.opts))
      }
      else {
        tmpf <- function(x) {
          objectivefunctiongr(x * psc)
        }
        oout$hessian <- try(jacobian(tmpf, oout$par/psc,
                                     method.args = hessian.opts))/outer(psc, psc)
      }
    }
  }
  if (skip.hessian || inherits(oout$hessian, "try-error"))
    oout$hessian <- namatrix
  coef <- oout$par
  nc <- names(coef)
  if (skip.hessian) {
    tvcov <- matrix(NA, length(coef), length(coef))
  }
  else {
    if (length(coef)) {
      if (use.ginv) {
        tmphess <- try(MASS::ginv(oout$hessian), silent = TRUE)
      }
      else {
        tmphess <- try(solve(oout$hessian, silent = TRUE))
      }
      if (inherits(tmphess, "try-error")) {
        tvcov <- matrix(NA, length(coef), length(coef))
        warning("couldn't invert Hessian")
      }
      else tvcov <- tmphess
    }
    else {
      tvcov <- matrix(numeric(0), 0, 0)
    }
  }
  dimnames(tvcov) <- list(nc, nc)
  min <- oout[[optimval]]
  fullcoef[nstart[order(oo)]] <- coef
  if (length(coef)) {
    gradvec <- if (!is.null(gr)) {
      objectivefunctiongr(coef)
    }
    else {
      if (inherits(tt <- try(grad(objectivefunction, coef),
                             silent = TRUE), "try-error"))
        NA
      else tt
    }
    oout$maxgrad <- max(abs(gradvec))
    if (!skip.hessian) {
      if (inherits(ev <- try(eigen(oout$hessian)$value,
                             silent = TRUE), "try-error"))
        ev <- NA
      oout$eratio <- min(Re(ev))/max(Re(ev))
    }
  }
  if (!is.null(conv <- oout$conv) && ((optimizer == "nlm" &&
                                       conv > 2) || (optimizer != "nlm" && conv != 0))) {
    if (is.null(oout$message)) {
      cmsg <- "unknown convergence failure: refer to optimizer documentation"
      if (optimizer == "optim") {
        if (conv == 1)
          cmsg <- "iteration limit 'maxit' reached"
        if (conv == 10)
          cmsg <- "degenerate Nelder-Mead simplex"
      }
      else if (optimizer == "nlm") {
        if (conv == 3)
          cmsg <- "last global step failed to locate a point lower than 'estimate': see ?nlm"
        if (conv == 4)
          cmsg <- "iteration limit exceeded"
        if (conv == 5)
          cmsg <- "maximum step size 'stepmax' exceeded five consecutive times: see ?nlm"
      }
    }
    else cmsg <- oout$message
    warning(paste0("convergence failure: code=", conv, " (",
                   cmsg, ")"))
  }
  m <- new("mle2", call = call, call.orig = call.orig, coef = coef,
           fullcoef = unlist(fullcoef), vcov = tvcov,
           details = oout, minuslogl = minuslogl, method = method,
           optimizer = optimizer, data = as.list(data), formula = formula)
  attr(m, "df") = length(m@coef)
  if (!missing(data))
    attr(m, "nobs") = length(data[[1]])
  environment(m) <- parent.frame()
  m
}
