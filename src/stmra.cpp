#include <Rcpp.h>
using namespace Rcpp;

// Compute a pairwise distance matrix of points in R2
// [[Rcpp::export]]
NumericMatrix R2_distances(NumericMatrix locA, NumericMatrix locB) {
  if (locA.ncol() != locB.ncol()) {
    // TODO: return error
  }
  NumericMatrix out(locA.nrow(), locB.nrow());
  for (uint32_t ia=0; ia<uint32_t(locA.nrow()); ++ia) {
    for (uint32_t ib=0; ib<uint32_t(locB.nrow()); ++ib) {
      double lonA = locA(ia,0);
      double latA = locA(ia,1) ;
      double lonB = locB(ib,0) ;
      double latB = locB(ib,1) ;
      out(ia, ib) = sqrt((lonA-lonB)*(lonA-lonB) + (latA-latB)*(latA-latB));
      if (std::isnan(out(ia, ib))) {
        out(ia, ib) = 0;
      }
    }
  }
  return out;
}

// Compute a pairwise distance matrix of points in R3
// [[Rcpp::export]]
NumericMatrix R3_distances(NumericMatrix locA, NumericMatrix locB) {
  if (locA.ncol() != locB.ncol()) {
    // TODO: return error
  }
  NumericMatrix out(locA.nrow(), locB.nrow());
  for (uint32_t ia=0; ia<uint32_t(locA.nrow()); ++ia) {
    for (uint32_t ib=0; ib<uint32_t(locB.nrow()); ++ib) {
      double xA = locA(ia,0);
      double yA = locA(ia,1) ;
      double zA = locA(ia,2) ;
      double xB = locB(ib,0);
      double yB = locB(ib,1) ;
      double zB = locB(ib,2) ;
      out(ia, ib) = sqrt((xA-xB)*(xA-xB)+(yA-yB)*(yA-yB)+(zA-zB)*(zA-zB));
      if (std::isnan(out(ia, ib))) {
        out(ia, ib) = 0;
      }
    }
  }
  return out;
}

// Compute a pairwise distance matrix of points on the surface of a sphere (S2) with radius 6371.0088 km
// [[Rcpp::export]]
NumericMatrix S2_distances(NumericMatrix locA, NumericMatrix locB) {
  if (locA.ncol() != locB.ncol()) {
    // TODO: return error
  }
  const double R = 6371.0088;
  NumericMatrix out(locA.nrow(), locB.nrow());
  for (uint32_t ia=0; ia<uint32_t(locA.nrow()); ++ia) {
    for (uint32_t ib=0; ib<uint32_t(locB.nrow()); ++ib) {
      double lonA = locA(ia,0) * PI / 180.0;
      double latA = locA(ia,1) * PI / 180.0;
      double lonB = locB(ib,0) * PI / 180.0;
      double latB = locB(ib,1) * PI / 180.0;
      out(ia, ib) = R*acos(sin(latA)*sin(latB) + cos(latA)*cos(latB)*cos(fabs(lonA-lonB)));
      if (std::isnan(out(ia, ib))) {
        out(ia, ib) = 0;
        //Rcpp::Rcout << ia << " " << ib << std::endl;
      }
    }
  }
  return out;
}


// Compute a pairwise distance matrix of points on the surface of a sphere (S2) as angles
// [[Rcpp::export]]
NumericMatrix S2_distances_angular(NumericMatrix locA, NumericMatrix locB) {
  if (locA.ncol() != locB.ncol()) {
    // TODO: return error
  }
  NumericMatrix out(locA.nrow(), locB.nrow());
  for (uint32_t ia=0; ia<uint32_t(locA.nrow()); ++ia) {
    for (uint32_t ib=0; ib<uint32_t(locB.nrow()); ++ib) {
      double lonA = locA(ia,0) * PI / 180.0;
      double latA = locA(ia,1) * PI / 180.0;
      double lonB = locB(ib,0) * PI / 180.0;
      double latB = locB(ib,1) * PI / 180.0;
      out(ia, ib) = acos(sin(latA)*sin(latB) + cos(latA)*cos(latB)*cos(fabs(lonA-lonB)));
      if (std::isnan(out(ia, ib))) {
        out(ia, ib) = 0;
        //Rcpp::Rcout << ia << " " << ib << std::endl;
      }
    }
  }
  return out;
}

// Compute a pairwise distance matrix of points in R1
// [[Rcpp::export]]
NumericMatrix R1_distances(NumericMatrix locA, NumericMatrix locB) {
  if (locA.ncol() != locB.ncol()) {
    // TODO: return error
  }
  NumericMatrix out(locA.nrow(), locB.nrow());
  for (uint32_t ia=0; ia<uint32_t(locA.nrow()); ++ia) {
    for (uint32_t ib=0; ib<uint32_t(locB.nrow()); ++ib) {
      out(ia, ib) = fabs(locA(ia,0) - locB(ib,0));
    }
  }
  return out;
}






double radial_basis_gauss_R1(double x, NumericVector w, double width = double(1)/double(20), uint16_t n_knots = 9, double lim0=-110, double lim1 = 110) {

  // TODO: assert that length(w) == n_knots
  std::vector<double> knots;
  for (uint16_t i=0; i < n_knots; ++i) {
    knots.push_back(lim0 + i* (lim1-lim0)/(n_knots-1)     );
  }
  double sum = 0.0;
  for (uint16_t i=0; i < n_knots; ++i) {
    double r = x - knots[i];
      double bf = exp(-(pow(width*r,2)));
    sum += w[i] * bf;
  }
  return sum;
}

// [[Rcpp::export]]
NumericMatrix cov_nonstationary_latitude(NumericMatrix locA, NumericMatrix locB, NumericVector theta, uint16_t n_knots = 9) {
  if (theta.length() != 3*n_knots) {
    Rcpp::stop("Invalid length of parameter vector theta");
  }
  if (n_knots != 9) {
    Rcpp::stop("Currently, the number of knots is fixed to 9 only.");
  }
  NumericMatrix out(locA.nrow(), locB.nrow());
  for (uint32_t ia=0; ia<uint32_t(locA.nrow()); ++ia) {
    for (uint32_t ib=0; ib<uint32_t(locB.nrow()); ++ib) {
      double xA = locA(ia,0);
      double yA = locA(ia,1) ;
      double xB = locB(ib,0) ;
      double yB = locB(ib,1) ;

      NumericVector weights_sd(n_knots);
      NumericVector weights_scale_x(n_knots);
      NumericVector weights_scale_y(n_knots);

      // extract weights from theta
      for (uint16_t i=0; i<n_knots; ++i) {
        weights_sd[i] = theta[i+0*n_knots];
        weights_scale_x[i] = theta[i+1*n_knots];
        weights_scale_y[i] = theta[i+2*n_knots];
      }

      double sd_1 = radial_basis_gauss_R1(yA, weights_sd);
      double sd_2 = radial_basis_gauss_R1(yB, weights_sd);
      double scale_x_1 = radial_basis_gauss_R1(yA, weights_scale_x);
      double scale_x_2 = radial_basis_gauss_R1(yB, weights_scale_x);
      double scale_y_1 = radial_basis_gauss_R1(yA, weights_scale_y);
      double scale_y_2 = radial_basis_gauss_R1(yB, weights_scale_y);
      double scale_x = (scale_x_1 + scale_x_2) / 2.0;
      double scale_y = (scale_y_1 + scale_y_2) / 2.0;
      double det = scale_x * scale_y;

      double dx = xA - xB;
      double dy = yA - yB;

      double C = dx*dx * (scale_y / det) + dy*dy * (scale_x / det);
      C = sd_1 * sd_2 * (1/sqrt(det)) * exp(-sqrt(C));

      out(ia, ib) = C;
      if (std::isnan(out(ia, ib))) {
        out(ia, ib) = 0;
      }
    }
  }
  return out;
}












