#include "clustRviz.h"

// [[Rcpp::export]]
Rcpp::List CARPL1_NF_FRAC(const arma::colvec& x,
                          int n,
                          int p,
                          double lambda_init,
                          double t,
                          const arma::colvec& weights,
                          const arma::colvec& uinit,
                          const arma::colvec& vinit,
                          const Eigen::SparseMatrix<double>& premat,
                          const arma::umat& IndMat,
                          const arma::umat& EOneIndMat,
                          const arma::umat& ETwoIndMat,
                          double rho = 1,
                          int max_iter = 1e4,
                          int burn_in = 50,
                          bool verbose=false,
                          int keep=10){


  int cardE = EOneIndMat.n_rows;
  arma::colvec unew(n*p);
  arma::colvec uold(n*p);
  unew = uinit;

  arma::colvec vnew(p*cardE);
  arma::colvec vold(p*cardE);
  vnew = vinit;

  arma::colvec lamnew(p*cardE);
  arma::colvec lamold(p*cardE);
  lamnew = vnew;

  arma::colvec proxin(p*cardE);
  double lambda = lambda_init;

  arma::uword path_iter = 0;
  arma::uword iter = 0;
  arma::mat UPath(n*p,1);
  UPath.col(iter)= uinit;
  arma::mat VPath(p*cardE,1);
  VPath.col(iter)= vinit;
  arma::colvec lambda_path(1);
  lambda_path(iter) = lambda_init;

  arma::colvec vZeroIndsnew(cardE,arma::fill::zeros);
  arma::colvec vZeroIndsold(cardE);
  arma::ucolvec vIdx(p);
  arma::mat vZeroInds_Path(cardE,1,arma::fill::zeros);

  arma::colvec arma_sparse_solver_input(n*p);
  int nzeros_old = 0;
  int nzeros_new = 0;
  Rcpp::List ret;


  while( (iter < max_iter) & (nzeros_new < cardE) ){
    iter += 1;
    uold = unew;
    vold = vnew;
    lamold = lamnew;
    vZeroIndsold = vZeroIndsnew;
    nzeros_old = nzeros_new;

    // u update
    // unew = arma::spsolve(premat, (1/rho)*x + (1/rho)*DtMatOpv2(rho*vold - lamold,n,p,IndMat,EOneIndMat,ETwoIndMat));
    arma_sparse_solver_input = (1/rho)*x + (1/rho)*DtMatOpv2(rho*vold - lamold,n,p,IndMat,EOneIndMat,ETwoIndMat);
    unew = cv_sparse_solve(premat, arma_sparse_solver_input);
    // v update
    proxin = DMatOpv2(unew,p,IndMat,EOneIndMat,ETwoIndMat) + (1/rho)*lamold;
    vnew = ProxL1(proxin,p,(1/rho)*lambda,weights);

    // lambda update
    lamnew = lamold + rho*(DMatOpv2(unew,p,IndMat,EOneIndMat,ETwoIndMat)-vnew);

    for(int l = 0; l<cardE; l++){
      vIdx = IndMat.row(l).t();
      if(sum(vnew.elem(vIdx)) == 0){
        vZeroIndsnew(l) = 1;
      }
    }
    nzeros_new = sum(vZeroIndsnew);
    if( (nzeros_new!=nzeros_old) | (iter % keep == 0)) {
      path_iter = path_iter + 1;
      UPath.insert_cols(path_iter,1);
      UPath.col(path_iter) = unew;
      VPath.insert_cols(path_iter,1);
      VPath.col(path_iter) = vnew;
      lambda_path.insert_rows(path_iter,1);
      lambda_path(path_iter) = lambda;

      vZeroInds_Path.insert_cols(path_iter,1);
      vZeroInds_Path.col(path_iter) = vZeroIndsnew;
    }

    if(iter >= burn_in){
      lambda = lambda*t;
    }


  }

  ret["u.path"] = UPath;
  ret["v.path"] = VPath;
  ret["v.zero.inds"]= vZeroInds_Path;
  ret["lambda.path"] = lambda_path;

  return(ret);
}
