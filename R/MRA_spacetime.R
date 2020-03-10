
#' Convert a raster stack to a data matrix
#'
#' @note This functions assumes equidistant images, i.e. a constant temporal gap between consecutive layers
#' @param s A RasterStack or RasterBrick object where layer dimension corresponds to time
#' @return A data.frame with rows representing pixels and columns for coordinates ("x","y","t") and data values ("values")
#' @export
stmra_stack_to_matrix <- function(s) {
  z = NULL
  for (i in 1:raster::nlayers(s)) {
    z = rbind(z, cbind(coordinates(s), i, as.vector(raster::subset(s, i))))
  }
  colnames(z) <- c("x", "y", "t", "value")
  return(z)
}

#' Create and initialize a spatiotemporal multiresolution approximation model
#' @param data A RasterBrick or RasterStack object with data used for parameter estimation
#' @param r length-one integer, number of knots in regions above the last partitioning level
#' @param M length-one integer, number of partitioning levels
#' @param covariance A valid spatiotemporal covariance function
#' @param region_minsize length-3 vector with minimum region size in x, y, and t dimensions (in this order)
#' @return A stmra object as a list with elements  TODO
#' @export
stmra_create_model <- function(data, M, r, cov_fun, region_minsize = c(0,0,0)) {

  resolution = c(NA, NA, NA)
  if ("RasterBrick" %in% class(data) | "RasterStack" %in% class(data))
  {
    z = stmra_stack_to_matrix(data)
    z = z[which(!is.na(z[,"value"])),] # ignore NA values
    domain = c("xmin" = extent(data)@xmin -0.0001, "xmax" = extent(data)@xmax + 0.0001,
               "ymin" = extent(data)@ymin -0.0001, "ymax" = extent(data)@ymax + 0.0001,
               "tmin" = 1-0.0001, "tmax" = nlayers(data) + 0.0001)
    resolution = c(raster::res(data), 1)
  }
  else if (is.matrix(data)) {
    # TODO: check for existing columns x, y, t, value
    z = data
    z = z[which(!is.na(z[,"value"])),] # ignore NA values
    domain = c("xmin" = min(z[,"x"],na.rm = TRUE) -0.0001, "xmax" = min(z[,"x"],na.rm = TRUE) + 0.0001,
               "ymin" = min(z[,"y"],na.rm = TRUE) -0.0001, "ymax" = min(z[,"y"],na.rm = TRUE) + 0.0001,
               "tmin" = min(z[,"t"],na.rm = TRUE) -0.0001, "tmax" = min(z[,"t"],na.rm = TRUE) + 0.0001)
  }
  else {
    stop("Expected Raster or matrix as data argument")
  }
  part.est = partition_spacetime(M=M, r=r, domain = domain,z = z[,"value"],locs=z[,c("x","y","t")],
                                 region_minsize = region_minsize)


  return(list(data = data,
              part_func = partition_spacetime,
              M = M,
              r = r,
              J = part.est$J,
              cov_fun = cov_fun,
              theta = NULL,
              part.est = part.est,
              region_minsize = region_minsize,
              resolution = resolution))
}


#' Estimate covariance function parameters of a spatiotemporal multiresolution approximation
#' @param model A stmra model as returned from \code{\link{stmra_create_model}}
#' @param theta0 initial values
#' @param lower_bounds lower boundaries of parameter values
#' @param upper_bounds upper boundaries of parameter values
#' @param ineq_constraints extended matrix with linear inequality constraints (rightmost column is RHS vector of the system)
#' @param control list to adjust optimization parameters e.g. maximum number of iterations, or tolerance (see Details)
#' @param optim.method optimization method (1 = BOBYQA from nloptr package, 2 = LL-BFGS-B, 3 = constrOptim from constrOptim package), see details
#' @param trace logical; print current parameter values and objective function value after each iteration
#' @details
#' The supported options and ther names in the control argument vary with the used optimization method. Please check the specific method for details.
#' @return model with estimated covariance function parameters theta
#' @export
stmra_estimate <- function(model, theta0, lower_bounds = NULL, upper_bounds = NULL, trace = FALSE, control = list(),
                           optim.method = 1, ineq_constraints = NULL) {



  part = model$part.est

  env_stmra_estimate = environment()

  model$estimation_log = data.frame()

  f = function(theta) {
    val = MRA.fast(theta, cov.fun = model$cov_fun,indices = part$indices,knots = part$knots, data = part$data, J=part$J)
    if(trace) {
      cat(paste("(#", nrow(env_stmra_estimate$model$estimation_log) + 1, ") THETA=(", paste(round(theta, digits=4),collapse = ", "), ") -> ", val, "\n", sep=""))
      env_stmra_estimate$model$estimation_log = rbind(env_stmra_estimate$model$estimation_log, cbind(t(theta), val))
    }
    return(val)
  }

  if (optim.method == 1) {
    if(!requireNamespace("nloptr", quietly = TRUE)) {
      stop("package nloptr not found; please install first")
    }

    if (is.null(lower_bounds)) {
      stop("missing lower bounds")
    }
    if (is.null(upper_bounds)) {
      stop("missing upper bounds")
    }
    if (!is.null(ineq_constraints)) {
      warning("ignoring inequality constraints, use optim.method = 3 if needed.")
    }
    res = nloptr::bobyqa(theta0, f, lower=lower_bounds, upper=upper_bounds, control=control)
  }
  else if (optim.method == 2) {
    if(!requireNamespace("nloptr", quietly = TRUE)) {
      stop("package nloptr not found; please install first")
    }
    if (is.null(lower_bounds)) {
      stop("missing lower bounds")
    }
    if (is.null(upper_bounds)) {
      stop("missing upper bounds")
    }
    if (!is.null(ineq_constraints)) {
      warning("ignoring inequality constraints, use optim.method = 3 if needed.")
    }
    res = nloptr::lbfgs(theta0, f, lower=lower_bounds, upper=upper_bounds, control=control)
  }
  else if (optim.method == 3) {
    if(!requireNamespace("constrOptim", quietly = TRUE)) {
      stop("package constrOptim not found; please install first")
    }

    if (!is.null(ftol.abs))
      warning("ftol.abs is not supportted for optim.method == 3 and will be ignored")

    np = length(theta0)
    ui = matrix(nrow=0,ncol=np)
    ci = numeric(0)
    if (!is.null(ineq_constraints)) {
      nc = nrow(ineq_constraints)
      if (ncol(ineq_constraints) != np + 1) {
        stop("ineq_constraints must have ncol = npars + 1")
      }
      ui = ineq_constraints[,1:np]
      ci = ineq_constraints[,np + 1]
    }

    if (!is.null(lower_bounds)) {
      ui = rbind(ui, diag(1, np)) # lower bound is included
      ci = c(ci, lower_bounds)
    }
    if (!is.null(upper_bounds)) {
      ui = rbind(ui, diag(-1, np))
      ci = c(ci, -upper_bounds) # upper bound is NOT included
    }
    res = constrOptim(theta0, f, NULL, ui =ui, ci = ci, control=control)
  }
  else {
    stop("Invalid optimization method")
  }
  if (nrow( model$estimation_log) > 0) {
    colnames(model$estimation_log) <- c(paste("theta_", 1:length(theta0), sep=""), "fun")
  }
  if (trace) {
    print(res)
  }
  model$theta = res$par
  return(model)
}



#' Predict values at unobserved locations
#'
#' @param model A spatiotemporal multiresolution approximation object including covariance function parameters theta
#' @param data Either a RasterBrick or RasterStack object, or a data matrix with named(!) columns "x", "y", "t", "value"; see Details
#' @param mask A RasterBrick / RasterStack object with same dimensions as `data`, only pixels with non NA cells in the mask will be predicted; ignored if data is a matrix
#' @param pred_locs matrix containing prediction locations with columns "x", "y", "t"
#' @return Depending on the type of `pred_locs` either a list of two RasterBricks with predicted means and variances, or a matrix with columns "x", "y", "t", "pred.mean",  and "pred.var"
#' @details
#'
#' The `data` argument may be either a data matrix, or a RasterBrick / RasterStack object. In the former case, the argument `pred_locs` must be provided. In the latter case, if neither `mask` nor `pred_locs` is provided, predictions are computed for all NA pixels in  `data`.
#' If furthermore an additional mask is given, only NA pixels with non-NA values in the mask are predicted. Both, `mask` and NA values in `data` are ignored if `pred_locs`` is provided.
#'
#' The returned type is controlled by the type of provided prediction locations. If `pred_locs` is a RasterStack / RasterBrick or NULL, the function returns a list with two RasterBricks as elements, containing
#' predicted means and variances repsectively.  If `pred_locs` is a matrix, the returned value as a matrix with the same number of rows but two more columns "pred.mean", and "pred.var", representing mean and variance predictions.k,
#'
#' @export
stmra_predict <- function(model, data, mask = NULL, pred_locs = NULL) {

  if (is.null(model$theta)) {
    stop("Model has no parameter estimates; call stmra_estimate() first")
  }
  if (is.null(model$resolution)) model$resolution = rep(NA,3)


  if (!is.null(pred_locs)) {
    if (!is.null(mask)) {
      warning("Using pred_locs as prediction locations, ignoring provided mask")
    }

    if ("RasterBrick" %in% class(data) | "RasterStack" %in% class(data)) {
      data = stmra_stack_to_matrix(data)
      domain = c("xmin" = extent(data)@xmin -0.0001, "xmax" = extent(data)@xmax + 0.0001,
                 "ymin" = extent(data)@ymin -0.0001, "ymax" = extent(data)@ymax + 0.0001,
                 "tmin" = 1-0.0001, "tmax" = nlayers(data) + 0.0001)
    }
    else {
      domain = c("xmin" = min(z[,"x"],na.rm = TRUE) -0.0001, "xmax" = min(z[,"x"],na.rm = TRUE) + 0.0001,
                 "ymin" = min(z[,"y"],na.rm = TRUE) -0.0001, "ymax" = min(z[,"y"],na.rm = TRUE) + 0.0001,
                 "tmin" = min(z[,"t"],na.rm = TRUE) -0.0001, "tmax" = min(z[,"t"],na.rm = TRUE) + 0.0001)
    }
    if ("RasterBrick" %in% class(pred_locs) | "RasterStack" %in% class(pred_locs)) {
      p = stmra_stack_to_matrix(pred_locs)
      pred_indexes = which(is.na(p[,"value"]))
      part =  model$part_func(M=model$M, r=model$r,domain = domain,
                              z = data[,"value"], locs=data[,c("x","y","t")],
                              pred_locs = p[pred_indexes,c("x","y","t")],
                              region_minsize = model$region_minsize,
                              resolution = model$resolution)
    }
    else {
      part =  model$part_func(M=model$M, r=model$r,domain = domain,
                              z = data[,"value"], locs=data[,c("x","y","t")],
                              pred_locs = pred_locs,
                              region_minsize = model$region_minsize,
                              resolution = model$resolution)
    }


  }
  else {
    if (!("RasterBrick" %in% class(data) | "RasterStack" %in% class(data))) {
      stop("Prediction locations are missing")
    }
    if(!is.null(mask)) {
      stopifnot(dim(data) == dim(mask))
      stopifnot(coordinates(data) == coordinates(mask))
    }

    z = stmra_stack_to_matrix(data)
    if (!is.null(mask)) {
      m = stmra_stack_to_matrix(mask)
      pred_indexes  = which(is.na(z[,"value"]) & !is.na(m[,"value"]))
    }
    else {
      pred_indexes = which(is.na(z[,"value"]))
    }
    part =  model$part_func(M=model$M, r=model$r,domain = c("xmin" = extent(data)@xmin -0.0001, "xmax" = extent(data)@xmax + 0.0001,
                                                            "ymin" = extent(data)@ymin -0.0001, "ymax" = extent(data)@ymax + 0.0001,
                                                            "tmin" = 1-0.0001, "tmax" = nlayers(data) + 0.0001),
                            z = z[-pred_indexes,"value"], locs=z[-pred_indexes,c("x","y","t")],
                            pred_locs = z[pred_indexes,c("x","y","t")],
                            region_minsize = model$region_minsize,
                            resolution = model$resolution)

  }

  pred <- MRA.fast(model$theta, cov.fun = model$cov_fun, indices = part$indices,knots = part$knots, data = part$data, pred.locs = part$pred.locs, J=part$J)

  xpred = NULL
  for (i in 1:length(pred)){
    if (!is.null(pred[[i]])) {
      xpred = rbind(xpred, cbind(part$pred.locs[[i]], pred[[i]]))
    }
  }

  colnames(xpred) = c("x", "y", "t", "pred.mean", "pred.var")
  xpred=xpred[order(xpred[,"t"],  xpred[,"y"], xpred[,"x"]),]


  # check that prediction location mapping is correct
  # assertthat::assert_that(all(xpred[,1] == z[pred_indexes[as.integer(row.names(xpred))], "x"]))
  # assertthat::assert_that(all(xpred[,2] == z[pred_indexes[as.integer(row.names(xpred))], "y"]))
  # assertthat::assert_that(all(xpred[,3] == z[pred_indexes[as.integer(row.names(xpred))], "t"]))

  # if input was raster, return raster as output
  if (!is.null(pred_locs)) {
    if ("RasterBrick" %in% class(pred_locs) | "RasterStack" %in% class(pred_locs)) {
      z.pred = rep(NA, prod(dim(pred_locs)))
      #z.pred[-pred_indexes[as.integer(row.names(xpred))]] = NA
      z.pred[pred_indexes[as.integer(row.names(xpred))]] <- xpred[,"pred.mean"]
      dim(z.pred) = dim(pred_locs)[c(2,1,3)]
      z.pred = aperm(z.pred, c(2,1,3))
      z.pred = raster::brick(z.pred)
      extent(z.pred) = extent(pred_locs)
      crs(z.pred) = crs(pred_locs)

      z.var= rep(NA, prod(dim(pred_locs)))
      #z.var[-pred_indexes[as.integer(row.names(xpred))]] = NA
      z.var[pred_indexes[as.integer(row.names(xpred))]] <- xpred[,"pred.var"]
      dim(z.var) = dim(pred_locs)[c(2,1,3)]
      z.var = aperm(z.var, c(2,1,3))
      z.var = raster::brick(z.var)
      extent(z.var) = extent(pred_locs)
      crs(z.var) = crs(pred_locs)

      return(list(pred.mean = z.pred, pred.var =  z.var))
    }
    else {
      return(xpred)
    }
  }
  else {
    # data must be RasterStack / RasterBrick
    z.pred = rep(NA, prod(dim(data)))
    #z.pred[-pred_indexes[as.integer(row.names(xpred))]] = NA
    z.pred[pred_indexes[as.integer(row.names(xpred))]] <- xpred[,"pred.mean"]
    dim(z.pred) = dim(data)[c(2,1,3)]
    z.pred = aperm(z.pred, c(2,1,3))
    z.pred = raster::brick(z.pred)
    extent(z.pred) = extent(data)
    crs(z.pred) = crs(data)

    z.var= rep(NA, prod(dim(data)))
    #z.var[-pred_indexes[as.integer(row.names(xpred))]] = NA
    z.var[pred_indexes[as.integer(row.names(xpred))]] <- xpred[,"pred.var"]
    dim(z.var) = dim(data)[c(2,1,3)]
    z.var = aperm(z.var, c(2,1,3))
    z.var = raster::brick(z.var)
    extent(z.var) = extent(data)
    crs(z.var) = crs(data)

    return(list(pred.mean = z.pred, pred.var =  z.var))
  }
}


#' Compare predictions with true observations
#'
#' @export
stmra_assess <- function(pred, truth) {


  out = list()

  pred.mean = as.vector(pred$pred.mean)
  pred_plus2sd  = pred.mean +  2 * sqrt(as.vector(pred$pred.var))
  pred_minus2sd = pred.mean -  2 * sqrt(as.vector(pred$pred.var))
  in2sd = as.vector(truth) > pred_minus2sd &  as.vector(truth) < pred_plus2sd

  out$measures = c(
    "RMSE" =  sqrt(mean(as.vector(truth - pred$pred.mean)^2, na.rm = TRUE)),
    "MSE"  =  mean(as.vector(truth - pred$pred.mean)^2, na.rm = TRUE),
    "MAE"  =  (mean(abs(as.vector(truth - pred$pred.mean)), na.rm = TRUE)),
    "COR"  = cor(as.vector(truth), as.vector(pred$pred.mean), use="complete.obs"),
    "WITHIN2SD" = sum(in2sd, na.rm = TRUE) / sum(!is.na(in2sd)),
    "R2" = 1 - (mean(as.vector(truth - pred$pred.mean)^2, na.rm = TRUE))/(mean((as.vector(truth) - mean(as.vector(truth), na.rm = TRUE))^2, na.rm = TRUE))
  )

  out$xyplot = function(...) {
    plot(x = as.vector(truth), y = as.vector(pred$pred.mean), xlab = "truth", ylab = "predicted", ...)
  }

  out$distribution_plot <- function(...) {
    true_values = as.vector(truth)
    pred_values = as.vector(pred$pred.mean)
    nonna_indexes = which(!is.na(pred_values) & !is.na(true_values))
    true_values = true_values[nonna_indexes]
    pred_values = pred_values[nonna_indexes]

    index_order = order(true_values,decreasing = TRUE)
    true_values = true_values[index_order]
    pred_values = pred_values[index_order]
    plot(x=1:length(true_values), y = true_values, type="h", ...)
    points(x=1:length(pred_values), y = pred_values, pch = ".", col = "red")
  }

  out$residuals = as.vector(truth) - as.vector(pred$pred.mean)
  return (out)
}







