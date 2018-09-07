#include "clustRviz.h"

// Take a vector of length n * k and re-order it as
// x(0), x(k), x(2*k), x(n*k), x(1), x(k + 1), x(2*k + 1), etc.
Eigen::VectorXd restride(const Eigen::VectorXd& x,
                         Eigen::Index k){
  Eigen::VectorXd ret(x.size());

  Eigen::Index len = x.size() / k;

  if(len * k != x.size()){
    Rcpp::stop("k does not divide the number of elements of x!");
  }

  // TODO: Optimize this!
  for(Eigen::Index j = 0; j < k; j++){
    for(Eigen::Index i = 0; i < len; i++){
      ret(i + j * len) = x(j + i * k);
    }
  }

  return ret;
}

// TODO - Document me!
Eigen::VectorXd DMatOpv2(const Eigen::VectorXd& u,
                         int p,
                         const Eigen::MatrixXi& IndMat,
                         const Eigen::MatrixXi& EOneIndMat,
                         const Eigen::MatrixXi& ETwoIndMat){

  Eigen::VectorXd out(EOneIndMat.rows() * p);

  // TODO - Optimize me!
  for(Eigen::Index i = 0; i < EOneIndMat.rows(); i++){
    Eigen::VectorXi out_index   = IndMat.row(i);
    Eigen::VectorXi e_one_index = EOneIndMat.row(i);
    Eigen::VectorXi e_two_index = ETwoIndMat.row(i);

    for(Eigen::Index j = 0; j < p; j++){
      out(out_index(j)) = u(e_one_index(j)) - u(e_two_index(j));
    }
  }

  return out;
}

// TODO -- Document me!
Eigen::VectorXd DtMatOpv2(const Eigen::VectorXd& v,
                          int n,
                          int p,
                          const Eigen::MatrixXi& IndMat,
                          const Eigen::MatrixXi& EOneIndMat,
                          const Eigen::MatrixXi& ETwoIndMat){

  Eigen::VectorXd out = Eigen::VectorXd::Zero(n * p);
  Eigen::VectorXi v_index(p);
  Eigen::VectorXi e_one_index(p);
  Eigen::VectorXi e_two_index;

  // TODO - Optimize me!
  for(Eigen::Index i = 0; i < IndMat.rows(); i++){
    v_index     = IndMat.row(i);
    e_one_index = EOneIndMat.row(i);
    e_two_index = ETwoIndMat.row(i);

    for(Eigen::Index j = 0; j < p; j++){
      out(e_one_index(j)) += v(v_index(j));
      out(e_two_index(j)) -= v(v_index(j));
    }
  }

  return out;
}

bool is_nan(double x){
  return x != x;
}

int sgn(double x){
  int ret;
  if(x == 0){
    ret = 0;
  } else{
    ret = x/std::abs(x);
  }
  return(ret);
}


// TODO -- Document me!
Eigen::VectorXd ProxL2(const Eigen::VectorXd& delta,
                       int p,
                       const Eigen::VectorXd& scalars,
                       const Eigen::MatrixXi& IndMat){

  Eigen::Index num_edges = scalars.size();
  Eigen::VectorXd ret(num_edges * p);

  if(p * num_edges != delta.size()){
    Rcpp::stop("p * num_edges != delta in ProxL2");
  }

  if(num_edges != IndMat.rows()){
    Rcpp::stop("Penalty weights not of same length as number of groups");
  }

  if(IndMat.size()!= delta.size()){
    Rcpp::stop("Do not have a group assignment for each data point.");
  }


  // TODO - Optimize this!
  for(Eigen::Index i = 0; i < num_edges; i++){
    Eigen::VectorXi index_vec = IndMat.row(i);

    Eigen::VectorXd this_delta = extract(delta, index_vec);
    double scale_factor = 1 - scalars(i) / this_delta.norm();

    if(is_nan(scale_factor)){
      for(Eigen::Index j = 0; j < this_delta.size(); j ++){
        ret(index_vec(j)) = this_delta(j);
      }
    } else if(scale_factor < 0){
      for(Eigen::Index j = 0; j < this_delta.size(); j ++){
        ret(index_vec(j)) = 0;
      }
    } else {
      for(Eigen::Index j = 0; j < this_delta.size(); j++){
        ret(index_vec(j)) = scale_factor * this_delta(j);
      }
    }
  }

  return ret;
}

double soft_thresh(double x, double lambda){
  if(std::abs(x) < lambda){
    return 0;
  } else if (x > 0){
    return x - lambda;
  } else {
    return x + lambda;
  }
}

// TODO - Document me!
Eigen::VectorXd ProxL1(const Eigen::VectorXd& delta,
                       int p,
                       double lambda,
                       const Eigen::VectorXd& weights){

  Eigen::Index num_edges = weights.size();
  Eigen::VectorXd ret(num_edges * p);

  // TODO - Optimize this!
  for(Eigen::Index i = 0; i < num_edges; i++){
    double weight = weights(i);

    for(Eigen::Index j = 0; j < p; j++){
      Eigen::Index elem_index = i * p + j;
      ret(elem_index) = soft_thresh(delta(elem_index), weight * lambda);
    }
  }

  return ret;
}
