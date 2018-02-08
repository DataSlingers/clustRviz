#' Compute CARP solution path
#'
#' \code{CARP} returns a fast approximation to the Convex Clustering
#' solution path along with visualizations such as dendrograms and
#' cluster paths. Visualizations may be static, interactive, or both.
#'
#' \code{CARP} solves the following problem.
#' \eqn{\textrm{min}_{U} \|X - U \| + \lambda \sum_{i <j} w_{ij} \| U_{\dot,i} - U_{\dot, j} \| }
#'
#' @param X A n by p matrix with rows the observations and columns the variables
#' @param obs.labels a vector of length n containing observations (row) labels
#' @param var.labels a vector of length p containing variable (column) labels
#' @param X.center A logical. Should X be centered?
#' @param X.scale A logical. Should X be scaled?
#' @param rho A positive number for augmented lagrangian. Not advisable to change.
#' @param phi A positive numner used for scaling in RBF kernel
#' @param weights A vector of positive number of length choose(n,2).
#' Determines observation pairs fusions weight.
#' @param k an integer >= 1. The number of neighbors used to create sparse weights
#' @param ncores an integer >= 1. The number of cores to use.
#' @param max.iter an integer. The maximum number of CARP iterations.
#' @param verbose either 0,1,or 2. Higher numbers indicate more verbosity while
#' computing.
#' @param burn.in an integer. The number of initial iterations at a fixed
#' value of (small) lambda_k
#' @param alg.type Which CARP algorithm to perform. Choices are 'carpviz'
#' and 'carp;
#' @param t a number greater than 1. The size of the multiplicitive
#' regularization parameter update. Typical values are: 1.1, 1.05, 1.01, 1.005
#' @param interactive A logical. Should interactive paths and dendrograms be
#' returned?
#' @param static A logical. Should static paths and dendrograms be returned?
#' @param npcs A integer >= 2. The number of principal components to compute
#' for path visualization.
#' @return X the original data matrix
#' @return carp.dend if static==TRUE, a dendrogam representation of the clustering solution path
#' @return carp.cluster.path.vis the CARP solution path
#' @return n.obs the number of observations
#' @return p.vars the number of variables
#' @return phi A positive numner used for scaling in RBF kernel
#' @return k an integer >= 1. The number of neighbors used to create sparse weights
#' @return burn.in an integer. The number of initial iterations at a fixed
#' value of (small) lambda_k
#' @return alg.type Which CARP algorithm to perform. Choices are 'carpviz'
#' and 'carp;
#' @return X.center A logical. Should X be centered?
#' @return X.scale A logical. Should X be scaled?
#' @examples
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X
#' carp.fit <- CARP(
#'     X=presidential_speech$X,
#'     obs.labels=presidential_speech$labels)
CARP <- function(X,
                 obs.labels=NULL,
                 X.center=TRUE,
                 X.scale=FALSE,
                 rho=1,
                 phi=1e-3,
                 weights=NULL,
                 k=NULL,
                 ncores=1,
                 max.iter=1e6,
                 burn.in=50,
                 verbose=1,
                 alg.type='carpviz',
                 t = 1.05,
                 interactive=TRUE,
                 static=TRUE,
                 npcs=4){
  require(tidyverse)
  require(cvxclustr)
  if(is.logical(verbose)){
    verbose.basic = TRUE
    verbose.deep=FALSE
  }else if(verbose == 1) {
    verbose.basic = TRUE
    verbose.deep = FALSE
  }else if(verbose == 2){
    verbose.basic = TRUE
    verbose.deep = TRUE
  }else{
    verbose.basic=FALSE
    verbose.deep=FALSE
  }




  # get labels
  if(is.null(obs.labels)){
    n.labels <- rownames(X)
  } else{
    n.labels <- obs.labels
  }

  # center and scale
  X.orig <- X
  if(X.center|X.scale){
    X %>%
      scale(center=X.center,scale=X.scale) %>%
      t() -> X
  } else{
    X <- t(X)
  }
  n.obs <- ncol(X)
  p.vars <- nrow(X)

  # get weights
  if(is.null(weights)){
    weights <- cvxclustr::kernel_weights(X,phi)
    if(is.null(k)){
       k <- MinKNN(weights,n.obs)
    }
    weights <- cvxclustr::knn_weights(weights,k,n.obs)
  }


  if(verbose.basic) cat("Pre-computing weight-based edge sets\n")
  PreCompList <- suppressMessages(ConvexClusteringPreCompute(
    X=X,
    weights = weights,
    ncores = ncores,
    rho=rho
  ))
  cardE <- nrow(PreCompList$E)

  if(verbose.basic) cat('Computing CARP Path\n')
  switch(
    alg.type,
    carpviz={
      CARPL2_VIS_FRAC(x=X[TRUE],
            n= as.integer(n.obs),
            p = as.integer(p.vars),
            lambda_init = 1e-8,
            weights = weights[weights!=0],
            uinit = as.matrix(PreCompList$uinit),
            vinit = as.matrix(PreCompList$vinit),
            premat = PreCompList$PreMat,
            IndMat = PreCompList$ind.mat,
            EOneIndMat = PreCompList$E1.ind.mat,
            ETwoIndMat = PreCompList$E2.ind.mat,
            rho = rho,
            max_iter = as.integer(max.iter),
            burn_in = as.integer(burn.in),
            verbose = verbose.deep,
            try_tol=1e-5,
            ti=10,
            t_switch=1.01,
            keep=1) -> carp.sol.path
    },
    carp={
      cat('carp\n')
    }
  )

  if(verbose.basic) cat('Post-processing\n')
  ISP(
    sp.path = carp.sol.path$v.zero.inds %>% t(),
    v.path = carp.sol.path$v.path,
    u.path = carp.sol.path$u.path,
    lambda.path = carp.sol.path$lambda.path,
    cardE = cardE
  ) -> carp.cluster.path
  carp.cluster.path$sp.path.inter %>% duplicated(fromLast=FALSE) -> sp.path.dups
  adj.path <- CreateAdjacencyPath(PreCompList$E,sp.path = carp.cluster.path$sp.path.inter,n.obs)
  clust.graph.path <- CreateClusterGraphPath(adj.path)
  clust.path <- GetClustersPath(clust.graph.path)
  clust.path.dups <- duplicated(clust.path,fromLast = FALSE)

  carp.cluster.path[['sp.path.dups']] <- sp.path.dups
  carp.cluster.path[['adj.path']] <- adj.path
  carp.cluster.path[['clust.graph.path']] <- clust.graph.path
  carp.cluster.path[['clust.path']] <- clust.path
  carp.cluster.path[['clust.path.dups']] <- clust.path.dups

  if(static|interactive){
    carp.dend <- CreateDendrogram(carp.cluster.path,n.labels)
  } else{
    carp.dend <- NULL
  }
  if(interactive){
    X.pca <- prcomp(t(X),scale. = FALSE,center = FALSE)
    X.pca.rot <- X.pca$rotation[,1:npcs]
    lapply(1:length(carp.cluster.path$clust.path),function(iter){
      U <- t(matrix(carp.cluster.path$u.path.inter[,iter],ncol=n.obs))%*%X.pca.rot
      names(U) <- paste('PC',1:npcs)
      U %>%
        as.data.frame() %>%
        tbl_df() %>%
        mutate(
          Iter = iter,
          Obs = 1:n(),
          Cluster = carp.cluster.path$clust.path[[iter]]$membership,
          Lambda = carp.cluster.path$lambda.path.inter[iter],
          ObsLabel = n.labels
        )
    }) %>%
      do.call(rbind.data.frame,.) %>%
      tbl_df() %>%
      group_by(Iter) %>%
      mutate(
        NCluster = length(unique(Cluster))
      ) %>%
      ungroup() %>%
      mutate(
        LambdaPercent = Lambda / max(Lambda)
      ) -> carp.cluster.path.vis
  } else{
    carp.cluster.path.vis <- NULL
  }
  carp.fit <- list(
    X = X.orig,
    carp.dend = carp.dend,
    carp.cluster.path.vis = carp.cluster.path.vis,
    n.obs = n.obs,
    p.vars = p.vars,
    phi = phi,
    k = k,
    burn.in = burn.in,
    alg.type = alg.type,
    X.center=X.center,
    X.scale=X.scale,
    static=static,
    interactive=interactive
  )
  class(carp.fit) <- 'CARP'
  return(carp.fit)
}

#' Print CARP Summary
#'
#' Prints a brief descripton of a fitted carp object.
#'
#' Reports number of observations and variables of dataset, any preprocessing
#' done by the \code{CARP} function, regularization weight information,
#' the type of CARP algorithm performed, and the visualizations returned.
#' @param carp.fit a CARP object returned by \code{CARP}
#' @examples
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X
#' carp.fit <- CARP(
#'     X=presidential_speech$X,
#'     obs.labels=presidential_speech$labels)
#' print(carp.fit)
print.CARP <- function(carp.fit){
  preprocess.string <- c('center','scale')

  switch(
    carp.fit$alg.type,
    carpviz={
      alg.string = 'CARP VIZ'
    },
    carp={
      alg.string = 'CARP'
    })
  viz.string <- c('Static Dend', 'Static Path','Interactive Dend/Path')
  cat('CARP Fit Summary\n')
  cat('Number of Observations:', carp.fit$n.obs,'\n')
  cat('Number of Variables:', carp.fit$p.vars,'\n')
  cat('Pre-processing:',preprocess.string[c(carp.fit$X.center,carp.fit$X.scale)],'\n')
  cat('Weights: RBF Kernel, phi =',carp.fit$phi, 'k =',carp.fit$k,'\n')
  cat('Algorithm:',alg.string,'\n')
  cat('Visualizations:',viz.string[c(carp.fit$static,carp.fit$static,carp.fit$interactive)],'\n')

  cat('Raw Data:\n')
  carp.fit$X[1:min(5,carp.fit$n.obs),1:min(5,carp.fit$p.vars)]

}
