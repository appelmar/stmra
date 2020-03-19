
#' Find duplicate knots in a spatiotemporal recursive partition
#' @param part spatiotemporal partition as returned from \code{\link{stmra_partition}}
#' @return A data.frame with rows corresponding to knots that are duplicate and columns including knot coordinates, a region identifier `region` and a unique index `num` of the knot within the region
duplicate_knots <- function(part) {
  knots = data.frame()
  for (i in 1:length(part$knots)) {
    if (is.null(part$knots[[i]])) next
    nknots = nrow(part$knots[[i]])
    if (nknots == 0) next
    level = length(part$indices[[i]])
    knots = rbind(knots, data.frame(region = rep(i, nknots), num= (1:nknots), level = rep(level, nknots), x =  part$knots[[i]][,"x"], y = part$knots[[i]][,"y"], t = part$knots[[i]][,"t"]))
  }
  colnames(knots) <- c("region", "num","m", "x","y","t")
  rownums = which(duplicated(knots[,c("x","y","t")]))

  return(knots[rownums,]) # check with nrow(duplicate_knots) > 0
}





#' Create a spaciotemporal recursive partition of a given domain of interest, and assign
#' knots, data, and prediction locations to regions
#'
#' @param r length-one integer, number of knots in regions above the last partitioning level
#' @param M length-one integer, number of partitioning levels
#' @param domain numeric length-6 vector including named(!) elements "xmin", "xmax", "ymin", "ymax", "tmin", "tmax", or RasterStack / RasterBrick
#' @param resolution size of pixels in x, y, t dimensions (in this order) or RasterStack / RasterBrick object, used to apply subpixel shifts while placing knots at upper partitioning levels.
#' @param region_minsize length-3 vector with minimum region size in x, y, and t dimensions (in this order)
#' @export
stmra_partition <-function(M,  r, domain, region_minsize = c(0, 0, 0), resolution = c(NA, NA, NA)) {

  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("package raster required, please install it first")
  }

  if ("RasterBrick" %in% class(resolution) | "RasterStack" %in% class(resolution)) {
    resolution = c(raster::res(resolution), 1)
  }
  else if ("RasterBrick" %in% class(domain) | "RasterStack" %in% class(domain)) {
    if (all(is.na(resolution))) {
      resolution = c(raster::res(domain), 1)
    }
  }
  if ("RasterBrick" %in% class(domain) | "RasterStack" %in% class(domain)) {
    domain = c("xmin" = raster::extent(domain)@xmin -0.0001, "xmax" = raster::extent(domain)@xmax + 0.0001,
               "ymin" = raster::extent(domain)@ymin -0.0001, "ymax" = raster::extent(domain)@ymax + 0.0001,
               "tmin" = 1-0.0001, "tmax" = raster::nlayers(domain) + 0.0001)
  }
  # TODO: add gdalcubes
  # else if ("cube" %in% domain) {
  #
  # }

  stopifnot(length(region_minsize) == 3)
  stopifnot(length(resolution) == 3)


  min_dx = ifelse(!is.na(resolution[1]), resolution[1], 0.001)
  min_dy = ifelse(!is.na(resolution[2]), resolution[2], 0.001)
  min_dt = ifelse(!is.na(resolution[3]), resolution[3],0.001)

  split_x = rep(TRUE, M)
  split_y = rep(TRUE, M)
  split_t = rep(TRUE, M)
  for (m in 1:M) {
    if (((domain["xmax"] - domain["xmin"])/(2^m)) < region_minsize[1]) {
      split_x[m] = FALSE
    }
    if (((domain["ymax"] - domain["ymin"])/(2^m)) < region_minsize[2]) {
      split_y[m] = FALSE
    }
    if (((domain["tmax"] - domain["tmin"])/(2^m)) < region_minsize[3]) {
      split_t[m] = FALSE
    }
    if (!split_x[m] & !split_y[m] & !split_t[m]) {
      warning(paste("M has been reduced to ", m, " due to values of region_minsize", sep = ""))
      M = m
      break
    }
  }

  J_all = NULL
  indices=list(list())
  ind.list=list()
  if(M>0) {
    for(m in 1:M){
      J = 2^sum(c(split_x[m],split_y[m],split_t[m]))
      J_all = c(J_all, J)
      # case J == 1 is already handled above
      ind.list=c(ind.list,list(1:J))
      indices=c(indices,as.list(data.frame(t(expand.grid(rev(ind.list))[,seq(m,1,by=-1)]))))
    }
  }
  n.ind= length(indices)

  # create knot tree structure
  knots=vector("list",n.ind)
  data=vector("list",n.ind)
  bounds=matrix(ncol=6,nrow=n.ind)
  colnames(bounds) <- c("xmin", "xmax", "ymin", "ymax", "tmin", "tmax")

  for (i in 1:n.ind) {
    m = length(indices[[i]])

    # 1. derive bounds
    bounds[i,] = domain[c("xmin", "xmax", "ymin", "ymax", "tmin", "tmax")]
    if (m > 0) {
      for (im in 1:m) {

        nsplits = sum(c(split_x[im], split_y[im], split_t[im]))

        if (nsplits == 3) {
          # make t,y,x all in c(0,1)
          t = floor((indices[[i]][im]-1)/4)
          y = floor((indices[[i]][im]-1)/2) %% 2
          x = (indices[[i]][im]-1) %% 2

          bounds[i,"xmin"] = bounds[i,"xmin"] + x * (bounds[i,"xmax"] - bounds[i,"xmin"]) / 2
          bounds[i,"xmax"] = bounds[i,"xmax"] + (x-1) * (bounds[i,"xmax"] - bounds[i,"xmin"]) / 2
          bounds[i,"ymin"] = bounds[i,"ymin"] + y * (bounds[i,"ymax"] - bounds[i,"ymin"]) / 2
          bounds[i,"ymax"] = bounds[i,"ymax"] + (y-1) * (bounds[i,"ymax"] - bounds[i,"ymin"]) / 2
          bounds[i,"tmin"] = bounds[i,"tmin"] + t * (bounds[i,"tmax"] - bounds[i,"tmin"]) / 2
          bounds[i,"tmax"] = bounds[i,"tmax"] + (t-1) * (bounds[i,"tmax"] - bounds[i,"tmin"]) / 2
        }
        else if (nsplits == 2) {
          if (!split_t[im]) {
            y = floor((indices[[i]][im]-1)/2)
            x = (indices[[i]][im]-1) %% 2
            bounds[i,"xmin"] = bounds[i,"xmin"] + x * (bounds[i,"xmax"] - bounds[i,"xmin"]) / 2
            bounds[i,"xmax"] = bounds[i,"xmax"] + (x-1) * (bounds[i,"xmax"] - bounds[i,"xmin"]) / 2
            bounds[i,"ymin"] = bounds[i,"ymin"] + y * (bounds[i,"ymax"] - bounds[i,"ymin"]) / 2
            bounds[i,"ymax"] = bounds[i,"ymax"] + (y-1) * (bounds[i,"ymax"] - bounds[i,"ymin"]) / 2
          }
          else if (!split_y[im]) {
            t = floor((indices[[i]][im]-1)/2)
            x = (indices[[i]][im]-1) %% 2
            bounds[i,"xmin"] = bounds[i,"xmin"] + x * (bounds[i,"xmax"] - bounds[i,"xmin"]) / 2
            bounds[i,"xmax"] = bounds[i,"xmax"] + (x-1) * (bounds[i,"xmax"] - bounds[i,"xmin"]) / 2
            bounds[i,"tmin"] = bounds[i,"tmin"] + t * (bounds[i,"tmax"] - bounds[i,"tmin"]) / 2
            bounds[i,"tmax"] = bounds[i,"tmax"] + (t-1) * (bounds[i,"tmax"] - bounds[i,"tmin"]) / 2
          }
          else if (!split_x[im]){
            t = floor((indices[[i]][im]-1)/2)
            y = (indices[[i]][im]-1) %% 2
            bounds[i,"ymin"] = bounds[i,"ymin"] + y * (bounds[i,"ymax"] - bounds[i,"ymin"]) / 2
            bounds[i,"ymax"] = bounds[i,"ymax"] + (y-1) * (bounds[i,"ymax"] - bounds[i,"ymin"]) / 2
            bounds[i,"tmin"] = bounds[i,"tmin"] + t * (bounds[i,"tmax"] - bounds[i,"tmin"]) / 2
            bounds[i,"tmax"] = bounds[i,"tmax"] + (t-1) * (bounds[i,"tmax"] - bounds[i,"tmin"]) / 2
          }

        }
        else if (nsplits == 1) {
          if (split_t[im]) {
            t = (indices[[i]][im]-1) %% 2
            bounds[i,"tmin"] = bounds[i,"tmin"] + t * (bounds[i,"tmax"] - bounds[i,"tmin"]) / 2
            bounds[i,"tmax"] = bounds[i,"tmax"] + (t-1) * (bounds[i,"tmax"] - bounds[i,"tmin"]) / 2
          }
          else if (split_y[im]) {
            y = (indices[[i]][im]-1) %% 2
            bounds[i,"ymin"] = bounds[i,"ymin"] + y * (bounds[i,"ymax"] - bounds[i,"ymin"]) / 2
            bounds[i,"ymax"] = bounds[i,"ymax"] + (y-1) * (bounds[i,"ymax"] - bounds[i,"ymin"]) / 2
          }
          else if (split_x[im]) {
            x = (indices[[i]][im]-1) %% 2
            bounds[i,"xmin"] = bounds[i,"xmin"] + x * (bounds[i,"xmax"] - bounds[i,"xmin"]) / 2
            bounds[i,"xmax"] = bounds[i,"xmax"] + (x-1) * (bounds[i,"xmax"] - bounds[i,"xmin"]) / 2
          }
        }
        else {
          stop(paste("No splits at m=", im, sep=""))
        }
      }
    }

    # 2. place knots and add data for lowest level
    if(length(indices[[i]])<M) {

      approx_r = ceiling(r^(1/3))^3
      approx_r_per_dim = approx_r^(1/3)

      # TODO: apply pixel sub-pixel shift
      # offset_t = sign(runif(1,-1,1)) * runif(1, 0, (bounds[i,"tmax"] - bounds[i,"tmin"])/(2*(approx_r_per_dim + 2)))
      # offset_x = sign(runif(1,-1,1)) * runif(1, 0, (bounds[i,"xmax"] - bounds[i,"xmin"])/(2*(approx_r_per_dim + 2)))
      # offset_y = sign(runif(1,-1,1)) * runif(1, 0, (bounds[i,"ymax"] - bounds[i,"ymin"])/(2*(approx_r_per_dim + 2)))

      # make sure that knots are never placed outside of the region
      dt = min(min_dt, (bounds[i,"tmax"] - bounds[i,"tmin"])/(2*(approx_r_per_dim + 2)))
      dx = min(min_dx, (bounds[i,"xmax"] - bounds[i,"xmin"])/(2*(approx_r_per_dim + 2)))
      dy = min(min_dy, (bounds[i,"ymax"] - bounds[i,"ymin"])/(2*(approx_r_per_dim + 2)))
      offset_t = sign(stats::runif(1,-1,1)) * stats::runif(1, 0, dt / 2)
      offset_x = sign(stats::runif(1,-1,1)) * stats::runif(1, 0, dx / 2)
      offset_y = sign(stats::runif(1,-1,1)) * stats::runif(1, 0, dy / 2)


      # TODO: what if dimensions are too small?
      x.grid = as.matrix(expand.grid(x = seq(bounds[i,"xmin"], bounds[i,"xmax"], length.out = approx_r_per_dim + 2)[2:(approx_r_per_dim+1)] + offset_x ,
                                     y = seq(bounds[i,"ymin"], bounds[i,"ymax"], length.out = approx_r_per_dim + 2)[2:(approx_r_per_dim+1)] + offset_y,
                                     t = seq(bounds[i,"tmin"], bounds[i,"tmax"], length.out = approx_r_per_dim + 2)[2:(approx_r_per_dim+1)] + offset_t) )


      # x.grid = data.frame(x = runif(r, bounds[i,"xmin"], bounds[i,"xmax"]),
      #                    y = runif(r, bounds[i,"ymin"], bounds[i,"ymax"]),
      #                    t = runif(r, bounds[i,"tmin"], bounds[i,"tmax"]))
      colnames(x.grid) = c("x", "y", "t")
      # knots[[i]] = as.matrix(x.grid)

      # sample from all knots to make sure we have exactly r knots
      knots[[i]] = x.grid[sample(1:nrow(x.grid), r), ]
      #data[[i]]=NULL # this removes list entry??
      data[[i]]=numeric(0)
    }
    else {
      # ind.sub=which(locs[,1]>bounds[i,"xmin"] & locs[,1]<=bounds[i,"xmax"] &
      #                 locs[,2]>bounds[i,"ymin"] & locs[,2]<=bounds[i,"ymax"] &
      #                 locs[,3]>bounds[i,"tmin"] & locs[,3]<=bounds[i,"tmax"])


      knots[[i]]=matrix(nrow=0,ncol=3)
      data[[i]]=numeric(0)
    }

  }

  pred.locs = NA
  stopifnot(length(knots) == length(indices) &&  length(knots) == length(data))
  p = list(indices=indices,knots=knots,data=data,bounds=bounds,pred.locs=pred.locs, J=J_all, M=M, r=r)
  dupl = duplicate_knots(p)

  if (nrow(dupl) > 0) {
    for (i in 1:nrow(dupl)) {
      kn = p$knots[[dupl[i, "region"]]]
      p$knots[[dupl[i, "region"]]] =  p$knots[[dupl[i, "region"]]][-dupl[i, "num"],]
      if (!is.null(p$data[[dupl[i, "region"]]])) {
        p$data[[dupl[i, "region"]]] =  p$data[[dupl[i, "region"]]][-dupl[i, "num"]]
      }
    }
    warning(paste(nrow(dupl), " duplicate knots have been removed from partition", sep=""))
  }
  class(p) <- "stmra_partition"
  return(p)
}






#' Add observations to an existing partition
#'
#' @param p existing partition as created with stmra_partition
#' @param locs a matrix with data locations (points are rows, coordinates columns) with named(!) columns "x", "y", and "t"
#' @param z numeric vector with data values and length identical to `nrows(locs)`
.partition_add_data <- function(p, locs, z) {
  stopifnot("stmra_partition" %in% class(p))
  stopifnot(is.matrix(locs) && ncol(locs) == 3)
  stopifnot(is.vector(z) && nrow(locs) == length(z))

  J = p$J
  M = length(p$indices[[length(p$indices)]])
  if (M == 0) {
    iregion = 1
  }
  else {
    iregion = num.ind(rep(1, M), J)
  }

  # for all regions in lowest partitioning level,add locs and z as knots
  for (i in iregion:length(p$indices)) {
    ind.sub=which(locs[,1]>p$bounds[i,"xmin"] & locs[,1]<=p$bounds[i,"xmax"] &
                    locs[,2]>p$bounds[i,"ymin"] & locs[,2]<=p$bounds[i,"ymax"] &
                    locs[,3]>p$bounds[i,"tmin"] & locs[,3]<=p$bounds[i,"tmax"])

    p$knots[[i]]=locs[ind.sub,,drop=FALSE]
    p$data[[i]]=z[ind.sub]
  }


  dupl = duplicate_knots(p)

  if (nrow(dupl) > 0) {
    for (i in 1:nrow(dupl)) {
      kn = p$knots[[dupl[i, "region"]]]
      p$knots[[dupl[i, "region"]]] =  p$knots[[dupl[i, "region"]]][-dupl[i, "num"],]
      if (!is.null(p$data[[dupl[i, "region"]]])) {
        p$data[[dupl[i, "region"]]] =  p$data[[dupl[i, "region"]]][-dupl[i, "num"]]
      }
    }
    warning(paste(nrow(dupl), " duplicate knots have been removed from partition", sep=""))
  }

  return(p)
}


#' Add prediction locations to an existing partition
#'
#' @param p existing partition as created with stmra_partition
#' @param locs a matrix with prediction locations (points are rows, coordinates columns) with named(!) columns "x", "y", and "t"
.partition_add_pred <- function(p, locs) {
  stopifnot("stmra_partition" %in% class(p))
  stopifnot(is.matrix(locs) && ncol(locs) == 3)
  J = p$J
  M = length(p$indices[[length(p$indices)]])
  n.ind= length(p$indices)
  p$pred.locs = p$knots

  if (M == 0) {
    iregion = 1
  }
  else {
    iregion = num.ind(rep(1, M), J)
  }

  for (i in iregion:n.ind) {
      ind.sub=which(locs[,1]>p$bounds[i,"xmin"] & locs[,1]<=p$bounds[i,"xmax"] &
                      locs[,2]>p$bounds[i,"ymin"] & locs[,2]<=p$bounds[i,"ymax"] &
                      locs[,3]>p$bounds[i,"tmin"] & locs[,3]<=p$bounds[i,"tmax"])
      p$pred.locs[[i]]=locs[ind.sub,,drop=FALSE]
      rownames(p$pred.locs[[i]]) = ind.sub
  }
  return(p)
}

#' @export
print.stmra_partition <- function(x, ...) {
  M = length(x$indices[[length(x$indices)]])

  cat("A spacetime partition with:\n")
  cat(paste0("... M = ", x$M, " partitioning levels\n"))
  cat(paste0("... r = ", x$r, " basis functions in regions m < M\n"))
  cat(paste0("... J = (", paste0(x$J, collapse=","), ") splits by level\n"))

  cat("Knots per level:\n")
  # knots per level
  K = list()
  C = list()
  for (i in 0:x$M) {
    K = c(K,0)
    C = c(C, 0)
  }
  for (i in 1:length(x$indices)) {
    m = length(x$indices[[i]])
    K[[m+1]] = K[[m+1]] + nrow(x$knots[[i]])
    C[[m+1]] = C[[m+1]] + 1
  }
  knots_per_level = data.frame(M = numeric(), total = numeric(), avg = numeric())
  for (i in 0:M) {
    knots_per_level = rbind(knots_per_level,  data.frame(M = i, total = K[[i+1]], avg = K[[i+1]] /  C[[i+1]] ))
  }
  print(knots_per_level, row.names = FALSE)

  # Observations
  n = 0
  for (i in num.ind(rep(1, M), x$J):length(x$indices)) {
    n = n + sum(!is.na(x$data[[i]]))
  }
  cat(paste0("Observations: ", n, "\n"))

  # Prediction locations
  n = 0
  if (is.list(x$pred_locs)) {
    for (i in num.ind(rep(1, M), x$J):length(x$indices)) {
      n = n + length(x$pred_locs[[i]])
    }
  }
  cat(paste0("Prediction locations: ", n, "\n"))
}






