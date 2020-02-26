
# locs1 = matrix with rows = points, columns = x,y,t coordinate
# theta = c(range, sill, nugget, scale)
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
    x=fields::rdist(locs1,locs2)
    exp(-x/theta[1]) *  (theta[2] + (abs(x)<1e-8)* theta[3])
  }
}


# locs1 = matrix with rows = points, columns = x,y,t coordinate
# theta = c(s_nugget, t_nugget, t_range, t_sill??, sd_lat[9], range_s_eastwest_lat[9], range_s_southnorth_minor_lat[9]
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
    ds= R2_distances(locs1, locs2)
    dt = R1_distances(locs1[,3, drop = FALSE], locs2[,3, drop = FALSE])
    theta[1] * exp(-ds/theta[2]) *  ((1-theta[4]) + (abs(ds)<1e-8)* theta[4]) * exp(-dt/theta[3]) *  (1-theta[5]+ (abs(dt)<1e-8)* theta[5])
  }
}



# theta = c(sill, range_space, range_time, nugget)
# see Porcu et al. (2016), eq. (15)
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



radial_basis_gauss_R1 <- function(x, w, width = 1/20, n.knots = 9, lim=c(-110,110)) {
  sapply(x, function(x) {
    knots = seq(lim[1], lim[2], length.out = n.knots)
    r = (x-knots)
    bf = exp(-(((width)*(r))^2))
    sum(w * bf)})
}

# theta = c(nugget_spatial, nugget_temporal, range_temporal, sd_weights[9], scale_x_weights[9], scale_y_weights[9])
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





