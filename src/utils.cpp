#include "clustRviz.h"

// Take a vector of length n * k and re-order it as
// x(0), x(k), x(2*k), x(n*k), x(1), x(k + 1), x(2*k + 1), etc.
//' @useDynLib clustRviz
arma::vec restride(const arma::vec& x,
                   arma::uword k){
  arma::vec ret(x.n_elem);

  arma::uword len = x.n_elem / k;

  if(len * k != x.n_elem){
    Rcpp::stop("k does not divide the number of elements of x!");
  }

  arma::uvec in_ix = arma::regspace<arma::uvec>(0, k, x.n_elem - 1);

  for(arma::uword i = 0; i < k; i++){
    ret.subvec(len * i, len * (i + 1) - 1) = x(in_ix + i);
  }

  return ret;
}

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

double TwoNorm(const arma::colvec& x){
  return std::pow(sum(x % x), 0.5);
}

double TwoNorm(const Eigen::VectorXd x){
  return x.norm();
}

Eigen::VectorXd cv_sparse_solve(const Eigen::SparseMatrix<double>& A,
                                const Eigen::VectorXd& b){

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
  solver.compute(A);

  return solver.solve(b);
}

arma::vec cv_sparse_solve(const Eigen::SparseMatrix<double>& A,
                          const arma::vec& b){

  const Eigen::VectorXd eigen_b = Eigen::Map<const Eigen::VectorXd>(b.memptr(),
                                                                    b.n_rows,
                                                                    1);

  Eigen::VectorXd solution_eigen = cv_sparse_solve(A, eigen_b);

  arma::mat solution_arma = arma::mat(solution_eigen.data(),
                                      solution_eigen.rows(),
                                      1,
                                      false,
                                      false);

  return solution_arma.col(0);
}

arma::colvec DMatOpv2(const arma::colvec& u,
                      int p,
                      const arma::umat& IndMat,
                      const arma::umat& EOneIndMat,
                      const arma::umat& ETwoIndMat){

  int cardE = EOneIndMat.n_rows;
  arma::ucolvec retIdx(p);
  arma::ucolvec EOneIdx(p);
  arma::ucolvec ETwoIdx(p);
  arma::colvec ret(cardE*p);

  for(int l = 0; l < cardE; l++){
    retIdx = IndMat.row(l).t();
    EOneIdx = EOneIndMat.row(l).t();
    ETwoIdx = ETwoIndMat.row(l).t();

    ret.elem(retIdx) = u.elem(EOneIdx) - u.elem(ETwoIdx);
  }

  return(ret);
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

arma::colvec DtMatOpv2(const arma::colvec& v,
                       int n,
                       int p,
                       const arma::umat& IndMat,
                       const arma::umat& EOneIndMat,
                       const arma::umat& ETwoIndMat){

  arma::colvec out(n*p,arma::fill::zeros);
  arma::ucolvec vIdx(p);
  arma::ucolvec EOneIdx(p);
  arma::ucolvec ETwoIdx(p);
  int cardE = IndMat.n_rows;

  for(int l = 0; l < cardE; l++){
    vIdx = IndMat.row(l).t();
    EOneIdx = EOneIndMat.row(l).t();
    ETwoIdx = ETwoIndMat.row(l).t();

    out.elem(EOneIdx) = out.elem(EOneIdx) + v.elem(vIdx);
    out.elem(ETwoIdx) = out.elem(ETwoIdx) - v.elem(vIdx);
  }

  return(out);
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

arma::colvec ProxL2(const arma::colvec& delta,
                    int p,
                    const arma::colvec& scalars,
                    const arma::umat& IndMat){

  int cardE = IndMat.n_rows;
  arma::ucolvec retIdx(p);

  arma::colvec valvec(2);
  valvec(0) = 0;

  arma::colvec ret(p*cardE);

  for(int l = 0; l<cardE;l++){
    retIdx = IndMat.row(l).t();
    valvec(1) = 1 - scalars(l)/TwoNorm(delta.elem(retIdx));

    if(valvec.has_nan()){
      ret.elem(retIdx) = delta.elem(retIdx);
    } else{
      ret.elem(retIdx) = valvec.max()*delta.elem(retIdx);
    }

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
    double scale_factor = 1 - scalars(i) / TwoNorm(this_delta);

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

arma::colvec ProxL1(const arma::colvec& delta,
                    int p,
                    double lambda,
                    const arma::colvec& weights){

  int cardE = weights.n_rows;
  arma::colvec theta(cardE*p);
  arma::colvec ret(cardE*p);
  arma::colvec valvec(2);
  double ols;
  arma::uword widx;
  valvec(0) = 0;
  for(int l = 0; l < cardE*p; l++){
    widx = std::floor(l/p);
    ols = delta(l)*weights(widx);
    valvec(1) = std::abs(ols) - (lambda*std::pow(weights(widx),2));
    theta(l) = (sgn(ols))*valvec.max();
    ret(l) = theta(l)/weights(widx);
  }
  return(ret);
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
