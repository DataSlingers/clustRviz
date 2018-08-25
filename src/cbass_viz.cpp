#include "clustRviz.h"

// [[Rcpp::export]]
Rcpp::List CBASS_VIS(const arma::colvec& x,
                     int n,
                     int p,
                     double lambda_init,
                     const arma::colvec& weights_col,
                     const arma::colvec& weights_row,
                     const arma::colvec& uinit_row,
                     const arma::colvec& uinit_col,
                     const arma::colvec& vinit_row,
                     const arma::colvec& vinit_col,
                     const Eigen::SparseMatrix<double>& premat_row,
                     const Eigen::SparseMatrix<double>& premat_col,
                     const arma::umat& IndMat_row,
                     const arma::umat& IndMat_col,
                     const arma::umat& EOneIndMat_row,
                     const arma::umat& EOneIndMat_col,
                     const arma::umat& ETwoIndMat_row,
                     const arma::umat& ETwoIndMat_col,
                     double rho         = 1,
                     int max_iter       = 10000,
                     int burn_in        = 50,
                     bool verbose       = false,
                     bool verbose_inner = false,
                     double try_tol     = 1e-3,
                     int ti             = 15,
                     double t_switch    = 1.01,
                     bool l1            = false){

  double t = 1.1;
  int try_iter;

  int cardE_row = EOneIndMat_row.n_rows;
  int cardE_col = EOneIndMat_col.n_rows;

  arma::colvec unew = x;
  arma::colvec uold(n*p);
  arma::colvec ut(n*p);

  arma::colvec vnew_col = vinit_col;
  arma::colvec vold_col(p*cardE_col);
  arma::colvec vnew_row = vinit_row;
  arma::colvec vold_row(n*cardE_row);

  arma::colvec lamnew_col = vnew_col;
  arma::colvec lamold_col(p*cardE_col);
  arma::colvec lamnew_row = vnew_row;
  arma::colvec lamold_row(n*cardE_row);

  arma::colvec proxin_row(n*cardE_row);
  arma::colvec proxin_col(p*cardE_col);

  arma::colvec yt(n*p);
  arma::colvec y(n*p);
  arma::colvec pnew(n*p,arma::fill::zeros);
  arma::colvec pold(n*p);
  arma::colvec pt(n*p);
  arma::colvec qnew(n*p,arma::fill::zeros);
  arma::colvec qold(n*p);

  double lambda = lambda_init;
  double lambda_old = lambda;
  double lambda_try_old;
  double lambda_lower;
  double lambda_upper;
  double try_eps;

  int nzerosold_row = 0;
  int nzerosnew_row = 0;
  int nzerosold_col = 0;
  int nzerosnew_col = 0;

  arma::uword path_iter = 0;
  arma::uword iter = 0;
  arma::mat UPath(n*p,1);
  UPath.col(path_iter)= uinit_col;
  arma::mat VPath_row(n*cardE_row,1);
  VPath_row.col(path_iter)= vinit_row;
  arma::mat VPath_col(p*cardE_col,1);
  VPath_col.col(path_iter)= vinit_col;
  arma::colvec lambda_path(1);
  lambda_path(path_iter) = lambda_init;

  arma::colvec vZeroIndsnew_row(cardE_row,arma::fill::zeros);
  arma::colvec vZeroIndsold_row(cardE_row);
  arma::ucolvec vIdx_row(n);
  arma::colvec vZeroIndsnew_col(cardE_col,arma::fill::zeros);
  arma::colvec vZeroIndsold_col(cardE_col,arma::fill::zeros);
  arma::ucolvec vIdx_col(p);

  arma::mat vZeroIndsPath_row(cardE_row,1,arma::fill::zeros);
  arma::mat vZeroIndsPath_col(cardE_col,1,arma::fill::zeros);

  arma::colvec arma_sparse_solver_input_row(n*p);
  arma::colvec arma_sparse_solver_input_col(n*p);

  Rcpp::List ret;
  bool rep_iter;

  while( ((nzerosnew_row < cardE_row) | (nzerosnew_col < cardE_col)) & (iter < max_iter)){
    iter = iter + 1;
    uold = unew;
    vold_row = vnew_row;
    vold_col = vnew_col;
    lamold_row = lamnew_row;
    lamold_col = lamnew_col;

    pold = pnew;
    qold = qnew;

    vZeroIndsold_row = vZeroIndsnew_row;
    vZeroIndsold_col = vZeroIndsnew_col;
    nzerosold_row = nzerosnew_row;
    nzerosold_col = nzerosnew_col;

    pt = restride(pold, p);
    ut = restride(uold, p);
    rep_iter = true;
    try_iter = 0;
    lambda_upper = lambda;
    lambda_lower = lambda_old;
    lambda_try_old = lambda;
    while(rep_iter){
      try_iter = try_iter + 1;
      ////////////// solve row problem
      // u update
      // yt = arma::spsolve(premat_row, (1/rho)*(ut+pt) + (1/rho)*DtMatOpv2(rho*vold_row - lamold_row,p,n,IndMat_row,EOneIndMat_row,ETwoIndMat_row));
      arma_sparse_solver_input_row = (1/rho)*(ut+pt) + (1/rho)*DtMatOpv2(rho*vold_row - lamold_row,p,n,IndMat_row,EOneIndMat_row,ETwoIndMat_row);
      yt = cv_sparse_solve(premat_row, arma_sparse_solver_input_row);
      // v update
      proxin_row = DMatOpv2(yt,n,IndMat_row,EOneIndMat_row,ETwoIndMat_row) + (1/rho)*lamold_row;
      if(l1){
        vnew_row = ProxL1(proxin_row, n, (1/rho) * lambda, weights_row);
      } else {
        vnew_row = ProxL2(proxin_row, n, (1/rho) * weights_row * lambda, IndMat_row);
      }

      // lambda update
      lamnew_row = lamold_row + rho*(DMatOpv2(yt,n,IndMat_row,EOneIndMat_row,ETwoIndMat_row)-vnew_row);
      ////////// end solve row problem

      y = restride(yt, n);
      pnew = uold + pold - y;
      ////////////// Solve col problem
      // u update
      // unew = arma::spsolve(premat_col, (1/rho)*(y + qold) + (1/rho)*DtMatOpv2(rho*vold_col - lamold_col,n,p,IndMat_col,EOneIndMat_col,ETwoIndMat_col));

      arma_sparse_solver_input_col = (1/rho)*(y + qold) + (1/rho)*DtMatOpv2(rho*vold_col - lamold_col,n,p,IndMat_col,EOneIndMat_col,ETwoIndMat_col);
      unew = cv_sparse_solve(premat_col, arma_sparse_solver_input_col);
      // v update
      proxin_col = DMatOpv2(unew,p,IndMat_col,EOneIndMat_col,ETwoIndMat_col) + (1/rho)*lamold_col;
      if(l1){
        vnew_col = ProxL1(proxin_col, p, (1/rho) * lambda, weights_col);
      } else {
        vnew_col = ProxL2(proxin_col, p, (1/rho) * lambda * weights_col, IndMat_col);
      }

      // lambda update
      lamnew_col = lamold_col + rho*(DMatOpv2(unew,p,IndMat_col,EOneIndMat_col,ETwoIndMat_col)-vnew_col);
      ////////////// End Solve col problem
      qnew = y + qold - unew;

      for(int l = 0; l<cardE_row; l++){
        vIdx_row = IndMat_row.row(l).t();
        if(sum(vnew_row.elem(vIdx_row)) == 0){
          vZeroIndsnew_row(l) = 1;
        }
      }
      for(int l = 0; l<cardE_col; l++){
        vIdx_col = IndMat_col.row(l).t();
        if(sum(vnew_col.elem(vIdx_col)) == 0){
          vZeroIndsnew_col(l) = 1;
        }
      }
      nzerosnew_row = sum(vZeroIndsnew_row);
      nzerosnew_col = sum(vZeroIndsnew_col);
      if( (nzerosnew_row == nzerosold_row) & (nzerosnew_col == nzerosold_col) & (try_iter == 1)){
        // same sparsity on first try
        // go to next iteration
        rep_iter = false;
      } else if((nzerosnew_row > nzerosold_row + 1) | (nzerosnew_col > nzerosold_col + 1) ){
        // too much sparsity
        // redo iteration
        vZeroIndsnew_col = vZeroIndsold_col;
        vZeroIndsnew_row = vZeroIndsold_row;
        if(try_iter == 1){
          lambda = 0.5*(lambda_lower + lambda_upper);
        } else{
          lambda_upper = lambda;
          lambda = (0.5)*(lambda_lower + lambda_upper);
        }
      } else if( (nzerosnew_row == nzerosold_row) & (nzerosnew_col == nzerosold_col) ){
        // not enough sparsity
        // redo iteration
        vZeroIndsnew_col = vZeroIndsold_col;
        vZeroIndsnew_row = vZeroIndsold_row;
        lambda_lower = lambda;
        lambda = (0.5)*(lambda_lower + lambda_upper);
      } else{
        // just enough sparsity
        rep_iter = false;
      }
      try_eps = std::abs(lambda - lambda_try_old) / std::abs(lambda_try_old);
      if(try_iter > ti){
        rep_iter = false;
      }


    }
    if( (nzerosnew_col > 0) | (nzerosnew_row >0) ){
      t = t_switch;
    }

    if((nzerosnew_row!=nzerosold_row) | (nzerosnew_col!=nzerosold_col)){
      path_iter = path_iter + 1;
      UPath.insert_cols(path_iter,1);
      UPath.col(path_iter) = unew;
      VPath_row.insert_cols(path_iter,1);
      VPath_row.col(path_iter) = vnew_row;
      VPath_col.insert_cols(path_iter,1);
      VPath_col.col(path_iter) = vnew_col;
      lambda_path.insert_rows(path_iter,1);
      lambda_path(path_iter) = lambda;
      lambda_old = lambda;

      vZeroIndsPath_row.insert_cols(path_iter,1);
      vZeroIndsPath_row.col(path_iter) = vZeroIndsnew_row;
      vZeroIndsPath_col.insert_cols(path_iter,1);
      vZeroIndsPath_col.col(path_iter) = vZeroIndsnew_col;
    }

    lambda = lambda*t;

  }

  ret["u.path"] = UPath;
  ret["v.row.path"] = VPath_row;
  ret["v.col.path"] = VPath_col;
  ret["v.row.zero.inds"] = vZeroIndsPath_row;
  ret["v.col.zero.inds"] = vZeroIndsPath_col;
  ret["lambda.path"] = lambda_path;

  return(ret);
}
