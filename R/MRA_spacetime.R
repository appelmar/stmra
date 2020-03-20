
#' Convert a raster stack to a data matrix
#'
#' @note This functions assumes equidistant images, i.e. a constant temporal gap between consecutive layers
#' @param s A RasterStack or RasterBrick object where layer dimension corresponds to time
#' @return A data.frame with rows representing pixels and columns for coordinates ("x","y","t") and data values ("values")
#' @export
stmra_stack_to_matrix <- function(s) {
  z = NULL
  for (i in 1:raster::nlayers(s)) {
    z = rbind(z, cbind(sp::coordinates(s), i, as.vector(raster::subset(s, i))))
  }
  colnames(z) <- c("x", "y", "t", "value")
  return(z)
}






#' Estimate covariance function parameters of a spatiotemporal multiresolution approximation
#' @param data A RasterBrick or RasterStack object or a matrix with data used for parameter estimation
#' @param cov_fun A valid spatiotemporal covariance function
#' @param partition A stmra partition as returned from \code{\link{stmra_partition}}
#' @param theta0 initial values
#' @param lower_bounds lower boundaries of parameter values
#' @param upper_bounds upper boundaries of parameter values
#' @param ineq_constraints extended matrix with linear inequality constraints (rightmost column is RHS vector of the system)
#' @param control list to adjust optimization parameters e.g. maximum number of iterations, or tolerance (see Details)
#' @param optim.method optimization method (1 = nloptr::bobyqa, 2 = nloptr::lbfgs, 3 = stats::constrOptim), see Details
#' @param trace logical; print current parameter values and objective function value after each iteration
#' @param theta vector with known parameter values, if provided, skip parameter estimation
#' @details
#' The supported options and ther names in the control argument vary with the used optimization method. Please check the specific method for details.
#' @return model with estimated covariance function parameters theta
#' @export
stmra <- function(partition, cov_fun, data, theta0, lower_bounds = NULL, upper_bounds = NULL, trace = FALSE, control = list(),
                           optim.method = 1, ineq_constraints = NULL, theta = NULL) {

  stopifnot("stmra_partition" %in% class(partition))


  model = list()
  model$part = partition
  model$cov_fun = cov_fun


  if (!is.null(theta)) {
    model$theta = theta
  }
  else {
    if ("RasterBrick" %in% class(data) | "RasterStack" %in% class(data))
    {
      z = stmra_stack_to_matrix(data)
      z = z[which(!is.na(z[,"value"])),] # ignore NA values
    }
    else if (is.matrix(data)) {
      # TODO: check for existing columns x, y, t, value
      z = data
      z = z[which(!is.na(z[,"value"])),] # ignore NA values
    }
    else {
      stop("Expected Raster or matrix as data argument")
    }
    part = .partition_add_data(model$part, z[,c("x","y","t")], z[,"value"])

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
        lower_bounds = rep(1e-5, length(theta0))
        warning(paste0("missing lower bounds, using (", paste0(lower_bounds, collapse=","), ")"))
      }
      if (is.null(upper_bounds)) {
        upper_bounds = rep(1e5, length(theta0))
        warning(paste0("missing upper bounds, using (", paste0(upper_bounds, collapse=","), ")"))
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
        lower_bounds = rep(1e-5, length(theta0))
        warning(paste0("missing lower bounds, using (", paste0(lower_bounds, collapse=","), ")"))
      }
      if (is.null(upper_bounds)) {
        upper_bounds = rep(1e5, length(theta0))
        warning(paste0("missing upper bounds, using (", paste0(upper_bounds, collapse=","), ")"))
      }
      if (!is.null(ineq_constraints)) {
        warning("ignoring inequality constraints, use optim.method = 3 if needed.")
      }
      res = stats::optim(theta0, f, method = "L-BFGS-B", lower = lower_bounds, upper=upper_bounds, control=control)
    }
    else if (optim.method == 3) {
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
      res = stats::constrOptim(theta0, f, NULL, ui =ui, ci = ci, control=control)
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
  }
  class(model) <- "stmra"
  return(model)
}



#' @export
print.stmra <- function(x, ...) {
  cat(paste0("Partition: ", "M = ", x$part$M, " | r = ", x$part$r, " | J = (", paste0(x$part$J, collapse=","),")\n" ))
  cat(paste0("Parameters: (", paste0(format(x$theta) , collapse = ","), ") \n"))
  L = NA
  if (!is.null(x$estimation_log)) {
    L =  format(min(x$estimation_log$fun))
  }
  cat(paste0("MRA log-likelihood: ",L, " \n"))
  # model domain
  domain = matrix(x$part$bounds[1,], nrow=2, ncol=3)
  colnames(domain) <- c("x","y","t")
  rownames(domain) <- c("min","max")
  cat(paste0("Domain:\n"))
  print(domain)
}





#' Predict values at unobserved locations
#'
#' @param object Model of class "stmra"
#' @param data Either a RasterBrick or RasterStack object, or a data matrix with named(!) columns "x", "y", "t", "value"; see Details
#' @param mask A RasterBrick / RasterStack object with same dimensions as `data`, only pixels with non NA cells in the mask will be predicted; ignored if data is a matrix
#' @param pred_locs matrix containing prediction locations with columns "x", "y", "t"
#' @param ... further arguments, will be ignored
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
predict.stmra <- function(object, ..., data, mask = NULL, pred_locs = NULL) {

  if (is.null(object$theta)) {
    stop("Model has no parameter estimates; call stmra() first")
  }


  if (!is.null(pred_locs)) {
    if (!is.null(mask)) {
      warning("Using pred_locs as prediction locations, ignoring provided mask")
    }

    if ("RasterBrick" %in% class(data) | "RasterStack" %in% class(data)) {
      data = stmra_stack_to_matrix(data)
    }
    data = data[which(!is.na(data[,"value"])),] # ignore NA values
    if ("RasterBrick" %in% class(pred_locs) | "RasterStack" %in% class(pred_locs)) {
      p = stmra_stack_to_matrix(pred_locs)
      pred_indexes = which(is.na(p[,"value"]))
      part = .partition_add_data(object$part, data[,c("x","y","t")], data[,"value"])
      part = .partition_add_pred(part, p[pred_indexes,c("x","y","t")])

    }
    else {
      part = .partition_add_data(object$part, data[,c("x","y","t")], data[,"value"])
      part = .partition_add_pred(part, pred_locs)
    }


  }
  else {
    if (!("RasterBrick" %in% class(data) | "RasterStack" %in% class(data))) {
      stop("Prediction locations are missing")
    }
    if(!is.null(mask)) {
      stopifnot(dim(data) == dim(mask))
      stopifnot(sp::coordinates(data) == sp::coordinates(mask))
    }

    z = stmra_stack_to_matrix(data)
    if (!is.null(mask)) {
      m = stmra_stack_to_matrix(mask)
      pred_indexes  = which(is.na(z[,"value"]) & !is.na(m[,"value"]))
      data_indexes  = which(!is.na(z[,"value"]))
      part = .partition_add_data(object$part, z[data_indexes,c("x","y","t")], z[data_indexes,"value"])
      part = .partition_add_pred(part, z[pred_indexes,c("x","y","t")])
    }
    else {
      pred_indexes = which(is.na(z[,"value"]))
      part = .partition_add_data(object$part, z[-pred_indexes,c("x","y","t")], z[-pred_indexes,"value"])
      part = .partition_add_pred(part, z[pred_indexes,c("x","y","t")])
    }


  }

  pred <- MRA.fast(object$theta, cov.fun = object$cov_fun, indices = part$indices,knots = part$knots, data = part$data, pred.locs = part$pred.locs, J=part$J)

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
      raster::extent(z.pred) = raster::extent(pred_locs)
      raster::crs(z.pred) = raster::crs(pred_locs)

      z.var= rep(NA, prod(dim(pred_locs)))
      #z.var[-pred_indexes[as.integer(row.names(xpred))]] = NA
      z.var[pred_indexes[as.integer(row.names(xpred))]] <- xpred[,"pred.var"]
      dim(z.var) = dim(pred_locs)[c(2,1,3)]
      z.var = aperm(z.var, c(2,1,3))
      z.var = raster::brick(z.var)
      raster::extent(z.var) = raster::extent(pred_locs)
      raster::crs(z.var) = raster::crs(pred_locs)
      out = list(pred.mean = z.pred, pred.var =  z.var)
      class(out) <- "stmra_prediction"
      return(out)
    }
    else {
      class(xpred) <- "stmra_prediction"
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
    raster::extent(z.pred) = raster::extent(data)
    raster::crs(z.pred) = raster::crs(data)

    z.var= rep(NA, prod(dim(data)))
    #z.var[-pred_indexes[as.integer(row.names(xpred))]] = NA
    z.var[pred_indexes[as.integer(row.names(xpred))]] <- xpred[,"pred.var"]
    dim(z.var) = dim(data)[c(2,1,3)]
    z.var = aperm(z.var, c(2,1,3))
    z.var = raster::brick(z.var)
    raster::extent(z.var) = raster::extent(data)
    raster::crs(z.var) = raster::crs(data)
    out = list(pred.mean = z.pred, pred.var =  z.var)
    class(out) <- "stmra_prediction"
    return(out)
  }
}



#' Plot stmra predictions
#'
#' @param x Prediction as returned by [predict]
#' @param ... Further plot arguments passed to raster::plot
#' @param variance logical; if TRUE, plot predicted variances instead of means
#' @export
plot.stmra_prediction <- function(x, ..., variance = FALSE) {

  if (is.matrix(x)) {
    # not yet implemented
    stop("not yet implemented for non-raster predictions")
  }
  else {
    stopifnot(is.list(x))
    stopifnot(!is.null(x$pred.mean))
    stopifnot(!is.null(x$pred.var))
    stopifnot("RasterBrick" %in% class(x$pred.mean) | "RasterStack" %in% class(x$pred.var))
    if (variance) {
      raster::plot(x$pred.var, ...)
    }
    else {
      raster::plot(x$pred.mean, ...)
    }
  }
}








#' Compare predictions with true observations
#'
#' @param pred Predictions as returned from \code{\link{predict.stmra}}
#' @param truth True observations as RasterStack / RasterBrick
#' @return List with computed prediction scores and residuals
#' @export
stmra_assess <- function(pred, truth) {

  stopifnot("stmra_prediction" %in% class(pred))


  out = list()
  if (is.matrix(pred)) {
    # TODO: check for correct columns in pred
    stopifnot(is.matrix(truth))
    stop("stmra_assess only supports raster predictions at the moment")
  }
  else { # raster
    stopifnot(is.list(pred))
    stopifnot(!is.null(pred$pred.mean))
    stopifnot(!is.null(pred$pred.var))
    if ("RasterBrick" %in% class(pred$pred.mean) | "RasterStack" %in% class(pred$pred.var)) {
      stopifnot("RasterBrick" %in% class(truth) | "RasterStack" %in% class(truth))
    }

    pred.mean = as.vector(pred$pred.mean)
    pred_plus2sd  = pred.mean +  2 * sqrt(as.vector(pred$pred.var))
    pred_minus2sd = pred.mean -  2 * sqrt(as.vector(pred$pred.var))
    in2sd = as.vector(truth) > pred_minus2sd &  as.vector(truth) < pred_plus2sd

    out$measures = c(
      "RMSE" =  sqrt(mean(as.vector(truth - pred$pred.mean)^2, na.rm = TRUE)),
      "MSE"  =  mean(as.vector(truth - pred$pred.mean)^2, na.rm = TRUE),
      "MAE"  =  (mean(abs(as.vector(truth - pred$pred.mean)), na.rm = TRUE)),
      "COR"  = stats::cor(as.vector(truth), as.vector(pred$pred.mean), use="complete.obs"),
      "COV2SD" = sum(in2sd, na.rm = TRUE) / sum(!is.na(in2sd)),
      "R2" = 1 - (mean(as.vector(truth - pred$pred.mean)^2, na.rm = TRUE))/(mean((as.vector(truth) - mean(as.vector(truth), na.rm = TRUE))^2, na.rm = TRUE))
    )

    out$residuals = as.vector(truth) - as.vector(pred$pred.mean)

    # out$xyplot = function(...) {
    #   plot(x = as.vector(truth), y = as.vector(pred$pred.mean), xlab = "truth", ylab = "predicted", ...)
    # }
    #
    # out$distribution_plot <- function(...) {
    #   true_values = as.vector(truth)
    #   pred_values = as.vector(pred$pred.mean)
    #   nonna_indexes = which(!is.na(pred_values) & !is.na(true_values))
    #   true_values = true_values[nonna_indexes]
    #   pred_values = pred_values[nonna_indexes]
    #
    #   index_order = order(true_values,decreasing = TRUE)
    #   true_values = true_values[index_order]
    #   pred_values = pred_values[index_order]
    #   plot(x=1:length(true_values), y = true_values, type="h", ...)
    #   points(x=1:length(pred_values), y = pred_values, pch = ".", col = "red")
    # }
  }
  class(out) <- "stmra_assessment"
  return (out)
}




#' @export
print.stmra_assessment <- function(x, ...) {
  cat("Prediction scores:\n")
  print(x$measures)
  cat("Residuals:\n")
  summary(x$residuals)
}





