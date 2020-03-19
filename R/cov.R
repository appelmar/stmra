
#' Metric exponential spatiotemporal covariance function
#'
#' Compute pairwise covariances between two sets of locations using
#' a metric exponential spatiotemporal model.
#'
#' @param locs1 matrix with rows representing points and columns representing coordinates
#' @param locs2 matrix with rows representing points and columns representing coordinates
#' @param theta vector of parameters (see Details)
#' @details
#' Parameter vector `theta` includes range, sill, nugget, spatiotemporal anisotropy (in this order).
#' @export
stmra_cov_metric_exp <- function(locs1, locs2, theta) {
  if (!is.matrix(locs1)) {
    # single point?
    locs1 = matrix(locs1,nrow=1)
  }
  if (!is.matrix(locs2)) {
    # single point?
    locs2 = matrix(locs2,nrow=1)
  }
  if( (nrow(locs1)==0) || (nrow(locs2)==0) ){
    0
  }
  else {
    locs1[,3] = theta[4] * locs1[,3]
    locs2[,3] = theta[4] * locs2[,3]
    x = R3_distances(locs1, locs2)
    exp(-x/theta[1]) *  (theta[2] + (abs(x)<1e-8)* theta[3])
  }
}


#' Separable exp / exp spatiotemporal covariance function
#'
#' Compute pairwise covariances between two sets of locations using
#' a separable exponential spatiotemporal model with two exponential functions.
#'
#' @param locs1 matrix with rows representing points and columns representing coordinates
#' @param locs2 matrix with rows representing points and columns representing coordinates
#' @param theta vector of parameters (see Details)
#' @details
#' Parameter vector `theta` includes joint sill, spatial range, temporal range, spatial nugget, temporal nugget  (in this order).
#' @export
stmra_cov_separable_exp <- function(locs1, locs2, theta) {
  if (!is.matrix(locs1)) {
    # single point?
    locs1 = matrix(locs1,nrow=1)
  }
  if (!is.matrix(locs2)) {
    # single point?
    locs2 = matrix(locs2,nrow=1)
  }
  if( (nrow(locs1)==0) || (nrow(locs2)==0) ){
    return(0)
  }
  else {
    ds = R2_distances(locs1, locs2)
    dt = R1_distances(locs1[,3, drop = FALSE], locs2[,3, drop = FALSE])
    theta[1] * exp(-ds/theta[2]) *  ((1-theta[4]) + (abs(ds)<1e-8)* theta[4]) * exp(-dt/theta[3]) *  (1-theta[5]+ (abs(dt)<1e-8)* theta[5])
  }
}



#' Non-separable spatiotemporal covariance function on spherical distances
#'
#' Compute pairwise covariances between two sets of locations using
#' a non-separable spatiotemporal model with Great circle distances from Porcu et al. (2016), eq. (15).
#'
#' @param locs1 matrix with rows representing points and columns representing coordinates
#' @param locs2 matrix with rows representing points and columns representing coordinates
#' @param theta vector of parameters (see Details)
#' @details
#' Parameter vector `theta` includes partial sill, spatial scale, temporal scale, nugget (in this order).
#' @references
#' Porcu, E., Bevilacqua, M., & Genton, M. G. (2016). Spatio-temporal covariance and cross-covariance functions of the great circle distance on a sphere. Journal of the American Statistical Association, 111(514), 888-898.
#' @export
stmra_cov_porcu_etal_15 <- function(locs1, locs2, theta) {
  if (!is.matrix(locs1)) {
    # single point?
    locs1 = matrix(locs1,nrow=1)
  }
  if (!is.matrix(locs2)) {
    # single point?
    locs2 = matrix(locs2,nrow=1)
  }
  if( (nrow(locs1)==0) || (nrow(locs2)==0) ){
    0
  }
  else {
    ds = S2_distances_angular(locs1,locs2)
    dt = R1_distances(locs1[,3, drop = FALSE], locs2[,3, drop = FALSE])
    #return(list(ds = ds, dt=dt))
    ((theta[1] + (dt < 1e6 & ds < 1e6)*theta[4]) / (1 + (6371.0088 * ds/ theta[2]))) * exp(-(dt)/(theta[3]*(1 + (6371.0088 * ds/ theta[2]))^(1/4)))
  }
}


# One-dimensional Gaussian radial basis function
radial_basis_gauss_R1 <- function(x, w, width = 1/20, n.knots = 9, lim=c(-110,110)) {
  sapply(x, function(x) {
    knots = seq(lim[1], lim[2], length.out = n.knots)
    r = (x-knots)
    bf = exp(-(((width)*(r))^2))
    sum(w * bf)})
}



#' Nonstationary spatiotemporal covariance function varying by latitude
#'
#' Compute pairwise covariances between two sets of locations using
#' a nonstationary kernel convolution approach with standard deviation and anisotropy varying as weighted sum of
#' Gaussian radial basis functions cenetered at different latitudes.
#'
#' @param locs1 matrix with rows representing points and columns representing coordinates
#' @param locs2 matrix with rows representing points and columns representing coordinates
#' @param theta vector of parameters (see Details)
#' @details
#' Parameter vector `theta` includes spatial nugget, temporal nugget, temporal range, standard deviation by latitude (9 values), spatial range in east-west direction (9 values), spatial range in south-north direction (9 values) (in this order).
#' @references
#' Paciorek, C. J., & Schervish, M. J. (2006). Spatial modelling using a new class of nonstationary covariance functions. Environmetrics: The official journal of the International Environmetrics Society, 17(5), 483-506.
#' @export
stmra_cov_nonstationary_R2_separable <- function(locs1, locs2, theta) {
  if (!is.matrix(locs1)) {
    # single point?
    locs1 = matrix(locs1,nrow=1)
  }
  if (!is.matrix(locs2)) {
    # single point?
    locs2 = matrix(locs2,nrow=1)
  }
  if( (nrow(locs1)==0) || (nrow(locs2)==0) ){
    0
  }
  else {
    ds = R2_distances(locs1, locs2)
    dt = R1_distances(locs1[,3, drop = FALSE], locs2[,3, drop = FALSE])
    Cspatial = cov_nonstationary_latitude(locs1, locs2, theta[4:length(theta)])
    C = ((1-theta[1]) + (abs(ds)<1e-8)* theta[1]) * Cspatial * exp(-dt/theta[3]) *  (1-theta[2]+ (abs(dt)<1e-8)* theta[2])
    #C =  Cspatial * exp(-dt/theta[3])
    return(C)
  }
}





