
#' Find duplicate knots in a spatiotemporal recursive partition
#' @param part spatiotemporal partition as returned from \code{\link{partition_spacetime}}
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
#' @param domain numeric length-6 vector including named(!) elements "xmin", "xmax", "ymin", "ymax", "tmin", "tmax"
#' @param locs a matrix with data locations (points are rows, coordinates columns) with named(!) columns "x", "y", and "t"
#' @param z numeric vector with data values and length identical to `nrows(locs)`
#' @param resolution size of pixels in x, y, t dimensions (in this order), used to apply subpixel shifts while placing knots at upper partitioning levels
#' @param pred_locs matrix of prediction locations with named(!) colums "x", "y", and "t", or `NULL` for parameter estimation
#' @param region_minsize length-3 vector with minimum region size in x, y, and t dimensions (in this order)
#' @export
partition_spacetime <-function(M,  r=ceiling(length(z)/8^M), domain, locs, z, pred_locs=NULL, region_minsize = c(0, 0, 0), resolution = c(NA, NA, NA)) {

  stopifnot(length(region_minsize) == 3)
  stopifnot(length(resolution) == 3)

  valid_indexes = which(!is.na(z))
  locs = locs[valid_indexes,,drop=FALSE]
  z = z[valid_indexes]

  min_dx = ifelse(!is.na(resolution[1]), resolution[1], min(diff(unique(sort(locs[,1]))),na.rm = TRUE))
  min_dy = ifelse(!is.na(resolution[2]), resolution[2], min(diff(unique(sort(locs[,2]))),na.rm = TRUE))
  min_dt = ifelse(!is.na(resolution[3]), resolution[3], min(diff(unique(sort(locs[,3]))),na.rm = TRUE))

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
      dx = min(min_dt, (bounds[i,"xmax"] - bounds[i,"xmin"])/(2*(approx_r_per_dim + 2)))
      dy = min(min_dt, (bounds[i,"ymax"] - bounds[i,"ymin"])/(2*(approx_r_per_dim + 2)))
      offset_t = sign(runif(1,-1,1)) * runif(1, 0, dt / 2)
      offset_x = sign(runif(1,-1,1)) * runif(1, 0, dx / 2)
      offset_y = sign(runif(1,-1,1)) * runif(1, 0, dy / 2)


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
    } else {
      ind.sub=which(locs[,1]>bounds[i,"xmin"] & locs[,1]<=bounds[i,"xmax"] &
                      locs[,2]>bounds[i,"ymin"] & locs[,2]<=bounds[i,"ymax"] &
                      locs[,3]>bounds[i,"tmin"] & locs[,3]<=bounds[i,"tmax"])


      knots[[i]]=locs[ind.sub,,drop=FALSE]
      data[[i]]=z[ind.sub]
    }

  }



  # prediction locations
  if(is.null(pred_locs)) pred.locs=NA  else {
    pred.locs=vector("list",n.ind)
    for(i in 1:n.ind) {
      if(length(indices[[i]])>=M) {

        ind.sub=which(pred_locs[,1]>bounds[i,"xmin"] & pred_locs[,1]<=bounds[i,"xmax"] &
                        pred_locs[,2]>bounds[i,"ymin"] & pred_locs[,2]<=bounds[i,"ymax"] &
                        pred_locs[,3]>bounds[i,"tmin"] & pred_locs[,3]<=bounds[i,"tmax"])
        pred.locs[[i]]=pred_locs[ind.sub,,drop=FALSE]
        rownames( pred.locs[[i]]) = ind.sub
      } else pred.locs[[i]]=knots[[i]]
    }
  }

  assertthat::assert_that(length(knots) == length(indices) &&  length(knots) == length(data))
  p = list(indices=indices,knots=knots,data=data,bounds=bounds,pred.locs=pred.locs, J=J_all)
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





