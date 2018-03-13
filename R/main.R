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
#' @importFrom dplyr %>%
#' @importFrom cvxclustr kernel_weights
#' @importFrom cvxclustr knn_weights
#' @importFrom dplyr n
#' @importFrom dplyr tbl_df
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @export
#' @examples
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X[1:10,1:4]
#' carp.fit <- CARP(
#'     X=Xdat,
#'     obs.labels=presidential_speech$labels[1:10])
CARP <- function(X,
                 obs.labels=NULL,
                 var.labels=NULL,
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
                 alg.type=c('carpviz','carp','carpl1','carpvizl1'),
                 t = 1.05,
                 interactive=TRUE,
                 static=TRUE,
                 npcs=4){
  alg.type <- match.arg(alg.type)
  Iter <- Cluster <- Lambda <- NULL
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

  if(is.null(var.labels)){
    p.labels <- colnames(X)
  } else{
    p.labels <- var.labels
  }
  colnames(X) <- p.labels
  rownames(X) <- n.labels

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
      CARPL2_NF_FRAC(x=X[TRUE],
                     n=as.integer(n.obs),
                     p = as.integer(p.vars),
                     lambda_init = 1e-8,
                     t = t,
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
                     verbose=verbose.deep,
                     keep=1) -> carp.sol.path
    },
    carpl1={
      CARPL1_NF_FRAC(x=X[TRUE],
                     n=as.integer(n.obs),
                     p = as.integer(p.vars),
                     lambda_init = 1e-8,
                     t = t,
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
                     verbose=verbose.deep,
                     keep=1) -> carp.sol.path
    },
    carpvizl1={
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
    X.pca <- stats::prcomp(t(X),scale. = FALSE,center = FALSE)
    X.pca.rot <- X.pca$rotation[,1:npcs]
    lapply(1:length(carp.cluster.path$clust.path),function(iter){
      U <- t(matrix(carp.cluster.path$u.path.inter[,iter],ncol=n.obs))%*%X.pca.rot
      names(U) <- paste('PC',1:npcs)
      U %>%
        as.data.frame() %>%
        dplyr::tbl_df() %>%
        dplyr::mutate(
          Iter = iter,
          Obs = 1:n(),
          Cluster = carp.cluster.path$clust.path[[iter]]$membership,
          Lambda = carp.cluster.path$lambda.path.inter[iter],
          ObsLabel = n.labels
        )
    }) %>%
      do.call(rbind.data.frame,.) %>%
      dplyr::tbl_df() %>%
      dplyr::group_by(Iter) %>%
      dplyr::mutate(
        NCluster = length(unique(Cluster))
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        LambdaPercent = Lambda / max(Lambda)
      ) -> carp.cluster.path.vis
  } else{
    carp.cluster.path.vis <- NULL
  }
  carp.fit <- list(
    X = X.orig,
    carp.dend = carp.dend,
    carp.cluster.path.vis = carp.cluster.path.vis,
    carp.sol.path = carp.sol.path,
    cardE = cardE,
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
#' @param x a CARP object returned by \code{CARP}
#' @param ... Unused additional generic arguements
#' @export
#' @examples
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X[1:10,1:4]
#' carp.fit <- CARP(
#'     X=Xdat,
#'     obs.labels=presidential_speech$labels[1:10])
#' print(carp.fit)
print.CARP <- function(x,...){
  preprocess.string <- c('center','scale')

  switch(
    x$alg.type,
    carpviz={
      alg.string = 'CARP-VIZ'
    },
    carp={
      alg.string = 'CARP'
    },
    carpl1={
      alg.string = 'CARP L1'
    },
    carpvizl1={
      alg.string = 'CARP-VIZ L1'
    }
  )
  viz.string <- c('Static Dend', 'Static Path','Interactive Dend/Path')
  cat('CARP Fit Summary\n')
  cat('Number of Observations:', x$n.obs,'\n')
  cat('Number of Variables:', x$p.vars,'\n')
  cat('Pre-processing:',preprocess.string[c(x$X.center,x$X.scale)],'\n')
  cat('Weights: RBF Kernel, phi =',x$phi, 'k =',x$k,'\n')
  cat('Algorithm:',alg.string,'\n')
  cat('Visualizations:',viz.string[c(x$static,x$static,x$interactive)],'\n')

  cat('Raw Data:\n')
  x$X[1:min(5,x$n.obs),1:min(5,x$p.vars)]

}

#' Plot function for CARP objects
#'
#' \code{plot.CARP} displays both static and interactive visualizations of
#' the CARP clustering solution path, including dendrograms and clustering
#' paths
#'
#' Possible visualizations of the CARP solution path are: (i) a static
#' dendrogram representing the clustering solution path at various
#' levels of regularization; (ii) a static clustering path representing
#' the coalescence of observations projected via PCA; and (iii) a
#' interactive dendorgram and clustering path visualization, showing
#' how the clustering solutions change with increased reguarliazation.
#' The first interactive tab shows a movie of real-time cluster solutions
#' while the second tabs shows the cluster solution for various numbers
#' of clusters
#' @param x a CARP object returned by \code{CARP}
#' @param type a string specifying the type of plot to produce. 'dendrogram'
#' produces the static cluster dendrogram; 'path' produces the static
#' cluster path; and 'interactive' produces an interactive visualization of
#' both the cluster path and dendrogram
#' @param axis a character vector of length two with elements as 'PC1','PC2',..
#' etc. Specifics which principal component axis to display for the 'path'
#' visualization
#' @param percent a number between 0 and 1. Specifies how far along the
#' CARP path the 'path' visualization should display.
#' @param max.nclust a positive integer. The maximum number of clusters
#' to display in the interactive plot.
#' @param min.nclust a positive value. The minimum number of clusters to
#' display in the interactive plot.
#' @param ... Unused additional generic arguements
#' @param blwd a positive number. Line width on dendrograms.
#' @param lcex a positive number. Label size on dendrograms.
#' @export
#' @importFrom shiny shinyApp
#' @importFrom shiny fluidPage
#' @importFrom shiny titlePanel
#' @importFrom shiny tabsetPanel
#' @importFrom shiny fluidRow
#' @importFrom shiny animationOptions
#' @importFrom shiny column
#' @importFrom shiny plotOutput
#' @importFrom shiny sliderInput
#' @importFrom shiny uiOutput
#' @importFrom shiny renderUI
#' @importFrom shiny tags
#' @importFrom shiny checkboxGroupInput
#' @importFrom shiny renderPlot
#' @importFrom stats as.dendrogram
#' @importFrom stats median
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr distinct
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr select_
#' @importFrom dplyr %>%
#' @importFrom grDevices adjustcolor
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X[1:10,1:4]
#' carp.fit <- CARP(
#'     X=Xdat,
#'     obs.labels=presidential_speech$labels[1:10])
#' plot(carp.fit,type='interactive')
#' }
plot.CARP <- function(
  x,
  type=c('dendrogram','path','interactive'),
  axis = c('PC1','PC2'),
  blwd=2,
  lcex=.6,
  percent=1,
  max.nclust=9,
  min.nclust=1,
  ...){

  type = match.arg(type)
  switch(
    type,
    dendrogram={
      x$carp.dend %>%
        stats::as.dendrogram() %>%
        dendextend::set("branches_lwd",blwd) %>%
        dendextend::set("labels_cex",lcex) %>%
        plot(ylab='Amount of Regularization')
    },
    path={
      plot.cols <- c(
        axis,
        'Iter',
        'Obs',
        'Cluster',
        'Lambda',
        'ObsLabel',
        'NCluster',
        'LambdaPercent'
      )
      plot.frame <- x$carp.cluster.path.vis[,plot.cols]
      names(plot.frame)[1:2] <- c('V1','V2')
      plot.frame %>%
        dplyr::filter(LambdaPercent <= percent) %>%
        dplyr::filter(Iter > x$burn.in) %>%
        ggplot2::ggplot(ggplot2::aes(x=V1,y=V2,group=Obs)) +
        ggplot2::geom_path(
          ggplot2::aes(x=V1,y=V2),
          linejoin = 'round',
          color='red',
          size=1
        )  +
        ggplot2::geom_point(
          ggplot2::aes(x=V1,y=V2),
          data=plot.frame %>% dplyr::filter(Iter==1),
          color='black',
          size=I(2)
        ) +
        ggrepel::geom_text_repel(
          ggplot2::aes(x=V1,y=V2,label=ObsLabel),
          size=I(3),
          data=plot.frame %>% dplyr::filter(Iter == 1)
        ) +
        ggplot2::guides(color=FALSE,size=FALSE) +
        ggplot2::theme(axis.title = ggplot2::element_text(size=15)) +
        ggplot2::theme(axis.text = ggplot2::element_text(size=10)) +
        ggplot2::xlab(axis[1]) +
        ggplot2::ylab(axis[2])

    },
    interactive={
      shiny::shinyApp(
        ui=shiny::fluidPage(
          shiny::tags$style(type="text/css",
                     ".recalculating { opacity: 1.0; }"
          ),
          shiny::titlePanel("Clustering Example"),
          shiny::tabsetPanel(
            shiny::tabPanel("Movie",
                     shiny::fluidRow(
                       shiny::column(2,
                              shiny::sliderInput(
                                "regcent_movie",
                                "Amount of Regularization",
                                min = 0,
                                max = 1,
                                value = 0,
                                # step=.07,animate = animationOptions(interval=700,loop=T)
                                step=.03,animate = shiny::animationOptions(interval=1000,loop=T)
                              ),
                              shiny::uiOutput("choose_columns")
                       ),
                       shiny::column(5,
                              shiny::plotOutput("pcapathplot_movie",height = "700px")#,width = '900px')
                       ),
                       shiny::column(5,
                              shiny::plotOutput("dendplot_movie",height='700px')
                       )
                     )
            ),
            shiny::tabPanel("Static",
                     shiny::fluidRow(
                       shiny::column(2,
                              shiny::sliderInput("regcent_static",
                                          "Number of Clusters",
                                          min = min.nclust,
                                          max = max.nclust,
                                          value = max.nclust,
                                          step=1),
                              shiny::uiOutput("choose_columns_static")

                       ),
                       shiny::column(5,
                              shiny::plotOutput("pcapathplot_static",height = "700px")#,width = '900px')
                       ),
                       shiny::column(5,
                              shiny::plotOutput("dendplot_static",height='700px')
                       )
                     )
            )
          )
        ),
        server=function(input,output){

          # Drop-down selection box for which data set
          # Check boxes
          output$choose_columns <- shiny::renderUI({

            # Get the data set with the appropriate name
            colnames <- paste('PC',1:4,sep='')

            # Create the checkboxes and select them all by default
            shiny::checkboxGroupInput("columns", "Choose columns",
                               choices  = colnames,
                               selected = colnames[1:2])
          })
          # Check boxes
          output$choose_columns_static <- shiny::renderUI({
            # Get the data set with the appropriate name
            colnames <- paste('PC',1:4,sep='')

            # Create the checkboxes and select them all by default
            shiny::checkboxGroupInput("columns_static", "Choose columns",
                               choices  = colnames,
                               selected = colnames[1:2])
          })


          output$dendplot_movie <- shiny::renderPlot({
            x$carp.cluster.path.vis %>%
              dplyr::filter(LambdaPercent <= input$regcent_movie)  %>%
              dplyr::select(NCluster) %>%
              unlist() %>%
              unname() %>%
              min -> ncl
            x$carp.dend %>%
              stats::as.dendrogram() %>%
              dendextend::set("branches_lwd",2) %>%
              dendextend::set("labels_cex",.6) %>%
              plot(ylab='Amount of Regularization',cex.lab=1.5)
            my.cols <- grDevices::adjustcolor(c('grey','black'),alpha.f = .2)
            my.rect.hclust(x$carp.dend,k=ncl,border=2,my.col.vec=my.cols,lwd=3)


          })

          output$pcapathplot_movie <- shiny::renderPlot({
            min.iter=5
            rename.list <- list(Obs = 'Obs',
                                Cluster = 'Cluster',
                                Iter = 'Iter',
                                ObsLabel = 'ObsLabel',
                                NCluster = 'NCluster',
                                Var1 = input$columns[1],
                                Var2 = input$columns[2])
            if(input$regcent_movie == 0){
              cl.iter = min.iter
            } else{
              x$carp.cluster.path.vis %>%
                dplyr::filter(LambdaPercent <= input$regcent_movie) %>%
                dplyr::select(Iter) %>%
                unlist() %>%
                unname()  %>%
                max() -> cl.iter
            }
            rename.list <- list(Obs = 'Obs',
                                Cluster = 'Cluster',
                                Iter = 'Iter',
                                ObsLabel = 'ObsLabel',
                                NCluster = 'NCluster',
                                Var1 = input$columns[1],
                                Var2 = input$columns[2])
            x$carp.cluster.path.vis %>%
              dplyr::select_(.dots=rename.list) -> x$carp.cluster.path.vis.rename

            x$carp.cluster.path.vis.rename %>%
              dplyr::filter(Iter == cl.iter) %>%
              dplyr::select(Obs,Cluster,Var1,Var2) %>%
              dplyr::rename(
                MaxVar1 = Var1,
                MaxVar2 = Var2
              ) -> cl.assgn
            x$carp.cluster.path.vis.rename %>%
              dplyr::filter(Iter <= cl.iter) %>%
              dplyr::left_join(
                cl.assgn,
                by=c('Obs')
              ) -> tmp
            tmp %>%
              dplyr::filter(Iter > min.iter ) %>%
              ggplot2::ggplot(aes(x=Var1,y=Var2,group=Obs)) +
              ggplot2::geom_path(
                ggplot2::aes(x=Var1,y=Var2),
                linejoin = 'round',
                color='red'
              ) +
              ggplot2::geom_point(
                ggplot2::aes(x=MaxVar1,y=MaxVar2),
                data = tmp %>% dplyr::filter(Iter==1),
                color='red',
                size=I(4)
              ) +
              ggplot2::geom_point(
                aes(x=Var1,y=Var2),
                data=tmp %>% dplyr::filter(Iter==1),
                color='black',
                size=I(4)
              ) +
              ggplot2::geom_text(
                ggplot2::aes(x=Var1,y=Var2,label=ObsLabel),
                size=I(6),
                data=tmp %>% dplyr::filter(Iter == 1)
              ) +
              ggplot2::guides(color=FALSE,size=FALSE) +
              ggplot2::theme(axis.title = ggplot2::element_text(size=25)) +
              ggplot2::theme(axis.text = ggplot2::element_text(size=20)) +
              ggplot2::xlab(input$columns[1]) +
              ggplot2::ylab(input$columns[2])

          })
          output$dendplot_static <- shiny::renderPlot({
            x$carp.dend %>%
              stats::as.dendrogram() %>%
              dendextend::set("branches_lwd",2) %>%
              dendextend::set("labels_cex",.6) %>%
              plot(ylab='Amount of Regularization',cex.lab=1.5)
            my.cols <- grDevices::adjustcolor(RColorBrewer::brewer.pal(n=input$regcent_static,'Set1'),alpha.f=.2)
            my.rect.hclust(x$carp.dend,k=input$regcent_static,border=2,my.col.vec=my.cols,lwd=3)
          })
          output$pcapathplot_static <- shiny::renderPlot({
            ncl <- input$regcent_static
            my.cols <- grDevices::adjustcolor(RColorBrewer::brewer.pal(n=ncl,'Set1'))[order(unique(stats::cutree(x$carp.dend,k=ncl)[x$carp.dend$order]))]
            x$carp.cluster.path.vis %>%
              dplyr::distinct(Iter,NCluster) %>%
              dplyr::filter(NCluster == ncl) %>%
              dplyr::select(Iter) %>%
              unlist() %>%
              unname() %>%
              stats::median() -> cl.iter
            cl.iter <- floor(cl.iter)
            rename.list <- list(Obs = 'Obs',
                                Cluster = 'Cluster',
                                Iter = 'Iter',
                                ObsLabel = 'ObsLabel',
                                NCluster = 'NCluster',
                                Var1 = input$columns_static[1],
                                Var2 = input$columns_static[2])
            x$carp.cluster.path.vis %>%
              dplyr::select_(.dots=rename.list) -> carp.cluster.path.vis.rename

            carp.cluster.path.vis.rename %>%
              dplyr::filter(Iter == cl.iter) %>%
              dplyr::select(Obs,Cluster,Var1,Var2) %>%
              dplyr::rename(
                MaxVar1 = Var1,
                MaxVar2 = Var2,
                PlotCluster=Cluster
              ) %>%
              dplyr::mutate(
                PlotCluster = as.factor(PlotCluster)
              ) -> cl.assgn

            carp.cluster.path.vis.rename %>%
              dplyr::filter(Iter <= cl.iter) %>%
              dplyr::left_join(
                cl.assgn,
                by=c('Obs')
              ) -> tmp
            tmp %>%
              dplyr::filter(Iter > 50 ) %>%
              ggplot2::ggplot(aes(x=Var1,y=Var2,group=Obs)) +
              ggplot2::geom_path(
                ggplot2::aes(x=Var1,y=Var2,color=PlotCluster),
                linejoin = 'round'
              ) +
              ggplot2::geom_point(
                ggplot2::aes(x=MaxVar1,y=MaxVar2,color=PlotCluster),
                data = tmp %>% dplyr::filter(Iter==1),
                size=I(4)
              ) +
              ggplot2::geom_point(
                ggplot2::aes(x=Var1,y=Var2),
                data=tmp %>% dplyr::filter(Iter==1),
                color='black',
                size=I(4)
              ) +
              ggplot2::geom_text(
                aes(x=Var1,y=Var2,label=ObsLabel),
                size=I(6),
                data=tmp %>% dplyr::filter(Iter == 1)
              ) +
              ggplot2::scale_color_manual(values=my.cols)+
              ggplot2::guides(color=FALSE,size=FALSE) +
              ggplot2::theme(axis.title = ggplot2::element_text(size=25)) +
              ggplot2::theme(axis.text = ggplot2::element_text(size=20)) +
              ggplot2::xlab(input$columns_static[1]) +
              ggplot2::ylab(input$columns_static[2])
          })



        }
        # End Shiny App
      )

    }
  )

}

#' Method for returning CARP and CBASS clustering solutions
#'
#' See \code{Clustering.CARP} and \code{Clustering.CBASS} for details
#'
#' @param x a CARP or CBASS object
#' @param ... additional arguements passed to \code{Clustering.CARP} or
#' \code{Clustering.CBASS}
#' @return CARP clustering solutions or CBASS biclustering solutions
#' @export
Clustering <- function(x,...) {
  UseMethod("Clustering", x)
}

#' Get clustering solution from a CARP object
#'
#' Returns cluster labels and cluster means at a point along the CARP path,
#' or the entire cluster sequence of cluster labels and means
#'
#' Passing either the desired number of clusters (\code{k}) or the percent
#' regularization (\code{percent}) returns the clustering assignment
#' and cluster means at the specific point along the CARP path. If neither
#' \code{k} nor \code{percent} are specified, all clusteirng assignments and
#' mean matricies are returned.
#'
#' @param x A CARP object returned by \code{CARP}
#' @param k An interger between 1 and \code{n.obs}. The number of unique
#' clusters
#' @param percent A number between 0 and 1. The percent of regularization at
#' which to cut the path.
#' @param ... Unused additional generic arguements
#' @return A list with elements
#' \describe{
#' \item{\code{clustering.assignment}}{
#' In the case where either \code{k} or \code{percent} is specified, a vector
#' of cluster labels of length \code{n.obs}.
#' In the case where neither \code{k} nor \code{percent} is specified, a
#' matrix of size \code{n.obs} by \code{n.obs}, each row specifying a unique
#' cluster assignment along the CARP path.
#' }
#' \item{\code{cluster.means}}{
#' In the case where either \code{k} or \code{percent} is specified, a matrix
#' of dimension \code{p.vars} by \code{length(unique(clustering.assignment))}
#' with each column a cluster mean.
#' In the case where neither \code{k} nor \code{percent} is specificed, a
#' list of matricies of length \code{n.obs}, representing the cluster means
#' for each cluster assignment along the CARP path.
#' }
#' }
#' @importFrom stats cutree
#' @export
#' @examples
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X[1:10,1:4]
#' carp.fit <- CARP(
#'     X=Xdat,
#'     obs.labels=presidential_speech$labels[1:10])
#' # Return the CARP iterate with k=5 clusters
#' carp.clustering <- Clustering(carp.fit,k=5)
#' # Examine the cluster labels
#' carp.clustering$clustering.assignment
#' # Examine the cluster means
#' head(carp.clustering$cluster.means)
#' # Return the whole sequence of solutions
#' carp.clustering.full <- Clustering(carp.fit)
#' # Examine the k=5 solution again
#' carp.clustering.full$clustering.assignment[5,]
#' # Examine the k=5 means again
#' head(carp.clustering.full$cluster.means[[5]])
Clustering.CARP <- function(x,k=NULL,percent=NULL,...){
  if(!is.null(k)){
      clust.assign <- stats::cutree(x$carp.dend,k=k)
      lapply(unique(clust.assign),function(cl.lab){
        apply(
          matrix(t(x$X)[,clust.assign==cl.lab],nrow=x$p.vars),
          1,
          mean
        )
      }) %>%
        do.call(cbind,.) -> clust.means
      clust.assign <- paste('cl',clust.assign,sep='')
      colnames(clust.means) <- unique(clust.assign)

  } else if(!is.null(percent)){
      clust.assign <- stats::cutree(x$carp.dend,h=percent)
      lapply(unique(clust.assign),function(cl.lab){
        apply(
          matrix(t(x$X)[,clust.assign==cl.lab],nrow=x$p.vars),
          1,
          mean
        )
      }) %>%
        do.call(cbind,.) -> clust.means
      clust.assign <- paste('cl',clust.assign,sep='')
      colnames(clust.means) <- unique(clust.assign)

  } else{
    lapply(1:x$n.obs,function(k){
      stats::cutree(x$carp.dend,k)
    }) %>%
      do.call(rbind,.) -> clust.assign
    apply(clust.assign,1,function(cl.ass){
      lapply(unique(cl.ass),function(cl.lab){
        apply(
          matrix(t(x$X)[,cl.ass==cl.lab],nrow=x$p.vars),
          1,
          mean
        )
      }) %>%
      do.call(cbind,.) -> tmp.means
      colnames(tmp.means) <- paste('cl',1:length(unique(cl.ass)),sep='')
      tmp.means
    }) -> clust.means
    clust.assign <- matrix(paste('cl',clust.assign,sep=''),nrow=x$n.obs)

  }
  list(
    clustering.assignment = clust.assign,
    cluster.means = clust.means
  )
}


#' Compute CBASS solution path
#'
#' \code{CBASS} returns a fast approximation to the Convex BiClustering
#' solution path along with visualizations such as dendrograms and
#' heatmaps. Visualizations may be static, interactive, or both.
#'
#' \code{CBASS} solves the following problem.
#'
#' @param X A n.obs by p.vars matrix with rows the observations and columns the variables
#' @param obs.labels a vector of length n.obs containing observations (row) labels
#' @param var.labels a vector of length p.vars containing variable (column) labels
#' @param X.center.global a logical. If TRUE, the global mean of X is removed.
#' @param rho A positive number for augmented lagrangian. Not advisable to change.
#' @param phi A positive numner used for scaling in RBF kernel
#' @param weights.obs A vector of positive number of length choose(n.obs,2).
#' Determines observation pair fusions weight.
#' @param weights.vars A vector of positive number of length choose(p.vars,2).
#' Determines variable pair fusions weight.
#' @param k.obs an integer >= 1. The number of neighbors used to create sparse
#' observation weights
#' @param k.var an integer >= 1. The number of neighbors used to create sparse
#' variable weights
#' @param ncores an integer >= 1. The number of cores to use.
#' @param max.iter an integer. The maximum number of CARP iterations.
#' @param verbose either 0,1,or 2. Higher numbers indicate more verbosity while
#' computing.
#' @param burn.in an integer. The number of initial iterations at a fixed
#' value of (small) lambda_k
#' @param alg.type Which CARP algorithm to perform. Choices are 'cbassviz'
#' and 'cbass';
#' @param t a number greater than 1. The size of the multiplicitive
#' regularization parameter update. Typical values are: 1.1, 1.05, 1.01, 1.005.
#' Not used CBASS-VIZ algorithms.
#' @param interactive A logical. Should an interactive heatmap be returned?
#' @param static A logical. Should observation and variable dendrograms
#' be returned?
#' @return X the original data matrix
#' @return cbass.cluster.path.obs the CBASS observation solution path
#' @return cbass.cluster.path.var the CBASS variable solution path
#' @return cbass.dend.obs if static==TRUE, a dendrogam representation of the
#' observation clustering solution path
#' @return cbass.dend.var if static==TRUE, a dendrogam representation of the
#' variable clustering solution path
#' @return phi.obs A positive numner used for scaling in observation RBF kernel
#' @return phi.var A positive numner used for scaling in variable RBF kernel
#' @return k.obs an integer >= 1. The number of neighbors used to create sparse
#' observation weights
#' @return k.var an integer >= 1. The number of neighbors used to create sparse
#' variable weights
#' @return burn.in an integer. The number of initial iterations at a fixed
#' value of (small) lambda_k
#' @return alg.type Which CARP algorithm to perform. Choices are 'cbassviz'
#' and 'cbass'
#' @return interactive A logical. Should an interactive heatmap be returned?
#' @return static A logical. Should observation and variable dendrograms
#' be returned?
#' @return obs.labels a vector of length n.obs containing observations (row) labels
#' @return var.labels a vector of length p.vars containing variable (column) labels
#' @return X.center.global a logical. If TRUE, the global mean of X is removed.
#' @importFrom cvxclustr knn_weights
#' @importFrom cvxclustr kernel_weights
#' @export
CBASS <- function(X,
                 obs.labels=NULL,
                 var.labels=NULL,
                 X.center.global=TRUE,
                 rho=1,
                 phi=1e-1,
                 weights.obs=NULL,
                 weights.vars=NULL,
                 k.obs=NULL,
                 k.var=NULL,
                 t=NULL,
                 ncores=1,
                 max.iter=1e6,
                 burn.in=50,
                 verbose=1,
                 alg.type='cbassviz',
                 interactive=TRUE,
                 static=TRUE){
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
  if(is.null(var.labels)){
    p.labels <- colnames(X)
  } else{
    p.labels <- var.labels
  }


  # center and scale
  colnames(X) <- var.labels
  rownames(X) <- obs.labels
  X.orig <- X
  if(X.center.global){
    X <- X - mean(X)
    X <- t(X)
  } else{
    X <- t(X)
  }
  n.obs <- ncol(X)
  p.vars <- nrow(X)


  phi.row=phi/n.obs
  weights.row <- cvxclustr::kernel_weights(t(X),phi.row)
  k.row <- MinKNN(weights.row,p.vars)
  weights.row <- cvxclustr::knn_weights(weights.row,k.row,p.vars)
  weights.row <- weights.row/sum(weights.row)
  weights.row <- weights.row/sqrt(n.obs)
  PreCompList.row <- suppressMessages(
    ConvexClusteringPreCompute(X=t(X),
                               weights = weights.row,
                               ncores = ncores,rho=rho)
  )
  cardE.row <- nrow(PreCompList.row$E)

  phi.col=phi/p.vars
  weights.col <- cvxclustr::kernel_weights(X,phi.col)
  k.col <- MinKNN(weights.col,n.obs)
  weights.col <- cvxclustr::knn_weights(weights.col,k.col,n.obs)
  weights.col <- weights.col/sum(weights.col)
  weights.col <- weights.col/sqrt(p.vars)
  PreCompList.col <- suppressMessages(
    ConvexClusteringPreCompute(X=X,
                               weights = weights.col,
                               ncores = ncores,rho=rho)
  )
  cardE.col <- nrow(PreCompList.col$E)



  if(verbose.basic) cat('Computing CBASS Path\n')
  switch(
    alg.type,
    cbassviz={
      BICARPL2_VIS(x=X[TRUE],
                   n= as.integer(n.obs),
                   p = as.integer(p.vars),
                   lambda_init = 1e-6,
                   weights_row = weights.row[weights.row!=0],
                   weights_col = weights.col[weights.col!=0],
                   uinit_row = as.matrix(PreCompList.row$uinit),
                   uinit_col = as.matrix(PreCompList.col$uinit),
                   vinit_row = as.matrix(PreCompList.row$vinit),
                   vinit_col = as.matrix(PreCompList.col$vinit),
                   premat_row = PreCompList.row$PreMat,
                   premat_col = PreCompList.col$PreMat,
                   IndMat_row = PreCompList.row$ind.mat,
                   IndMat_col = PreCompList.col$ind.mat,
                   EOneIndMat_row = PreCompList.row$E1.ind.mat,
                   EOneIndMat_col = PreCompList.col$E1.ind.mat,
                   ETwoIndMat_row = PreCompList.row$E2.ind.mat,
                   ETwoIndMat_col = PreCompList.col$E2.ind.mat,
                   rho = rho,
                   max_iter = as.integer(max.iter),
                   burn_in = burn.in,
                   verbose = verbose.deep,
                   verbose_inner = verbose.deep,
                   try_tol = 1e-3,
                   ti = 10,
                   t_switch = 1.01) -> bicarp.sol.path
    },
    cbass={
      cat('carp\n')
    }
  )

  if(verbose.basic) cat('Post-processing\n')
  ISP(
    sp.path = bicarp.sol.path$v.row.zero.inds %>% t(),
    v.path = bicarp.sol.path$v.row.path,
    u.path = bicarp.sol.path$u.path,
    lambda.path = bicarp.sol.path$lambda.path,
    cardE = sum(weights.row != 0)
  ) -> bicarp.cluster.path.row
  bicarp.cluster.path.row$sp.path.inter %>% duplicated(fromLast=FALSE) -> sp.path.dups.row
  adj.path.row <- CreateAdjacencyPath(PreCompList.row$E,sp.path = bicarp.cluster.path.row$sp.path.inter,p.vars)
  clust.graph.path.row <- CreateClusterGraphPath(adj.path.row)
  clust.path.row <- GetClustersPath(clust.graph.path.row)
  clust.path.dups.row <- duplicated(clust.path.row,fromLast = FALSE)

  bicarp.cluster.path.row[['sp.path.dups']] <- sp.path.dups.row
  bicarp.cluster.path.row[['adj.path']] <- adj.path.row
  bicarp.cluster.path.row[['clust.graph.path']] <- clust.graph.path.row
  bicarp.cluster.path.row[['clust.path']] <- clust.path.row
  bicarp.cluster.path.row[['clust.path.dups']] <- clust.path.dups.row

  bicarp.dend.row <- CreateDendrogram(bicarp.cluster.path.row,p.labels)



  ISP(
    sp.path = bicarp.sol.path$v.col.zero.inds %>% t(),
    v.path = bicarp.sol.path$v.col.path,
    u.path = bicarp.sol.path$u.path,
    lambda.path = bicarp.sol.path$lambda.path,
    cardE = sum(weights.col != 0)
  ) -> bicarp.cluster.path.col
  bicarp.cluster.path.col$sp.path.inter %>% duplicated(fromLast=FALSE) -> sp.path.dups.col
  adj.path.col <- CreateAdjacencyPath(PreCompList.col$E,sp.path = bicarp.cluster.path.col$sp.path.inter,n.obs)
  clust.graph.path.col <- CreateClusterGraphPath(adj.path.col)
  clust.path.col <- GetClustersPath(clust.graph.path.col)
  clust.path.dups.col <- duplicated(clust.path.col,fromLast = FALSE)

  bicarp.cluster.path.col[['sp.path.dups']] <- sp.path.dups.col
  bicarp.cluster.path.col[['adj.path']] <- adj.path.col
  bicarp.cluster.path.col[['clust.graph.path']] <- clust.graph.path.col
  bicarp.cluster.path.col[['clust.path']] <- clust.path.col
  bicarp.cluster.path.col[['clust.path.dups']] <- clust.path.dups.col

  bicarp.dend.col <- CreateDendrogram(bicarp.cluster.path.col,n.labels)
  cbass.fit <- list(
    X = X.orig,
    cbass.sol.path = bicarp.sol.path,
    cbass.cluster.path.obs = bicarp.cluster.path.col,
    cbass.cluster.path.var = bicarp.cluster.path.row,
    cbass.dend.var = bicarp.dend.row,
    cbass.dend.obs = bicarp.dend.col,
    n.obs = n.obs,
    p.vars = p.vars,
    phi.var = phi.row,
    phi.obs = phi.col,
    k.obs = k.col,
    k.var = k.row,
    burn.in = burn.in,
    alg.type = alg.type,
    X.center.global=X.center.global,
    static=static,
    interactive=interactive,
    obs.labels = n.labels,
    var.labels = p.labels
  )
  class(cbass.fit) <- 'CBASS'
  return(cbass.fit)

}

#' Print CBASS Summary
#'
#' Prints a brief descripton of a fitted cbass object.
#'
#' Reports number of observations and variables of dataset, any preprocessing
#' done by the \code{CBASS} function, regularization weight information,
#' the type of CBASS algorithm performed, and the visualizations returned.
#' @param x a CBASS object returned by \code{CARP}
#' @param ... Unused additional generic arguements
#' @export
#' @examples
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X[1:10,1:4]
#' cbass.fit <- CBASS(
#'     X=Xdat)
#' print(cbass.fit)
print.CBASS <- function(x,...){
  preprocess.string <- c('global-center')

  switch(
    x$alg.type,
    cbassviz={
      alg.string = 'CBASS VIZ'
    },
    cbass={
      alg.string = 'CBASS'
    })
  viz.string <- c('Static Dend', 'Static Heatmap','Interactive Heatmap')
  cat('CBASS Fit Summary\n')
  cat('Number of Observations:', x$n.obs,'\n')
  cat('Number of Variables:', x$p.vars,'\n')
  cat('Pre-processing:',preprocess.string[c(x$X.center.global)],'\n')
  cat('Obs. Weights: RBF Kernel, phi =',x$phi.obs, 'k =',x$k.obs,'\n')
  cat('Var. Weights: RBF Kernel, phi =',x$phi.var, 'k =',x$k.var,'\n')
  cat('Algorithm:',alg.string,'\n')
  cat('Visualizations:',viz.string[c(x$static,x$static,x$interactive)],'\n')

  cat('Raw Data:\n')
  x$X[1:min(5,x$n.obs),1:min(5,x$p.vars)]

}


#' Plot function for CBASS objects
#'
#' \code{plot.CBASS} displays both static and interactive visualizations of
#' the CBASS clustering solution path, including dendrograms and heatmaps.
#'
#' Possible visualizations of the CBASS solution path are: (i) a static
#' dendrogram representing the clustering solution path at various
#' levels of regularization (both observations and variables); and (ii) an
#' interactive cluster heatmap, showing
#' how the clustering solutions change with increased reguarliazation.
#' @param x a CBASS object returned by \code{CARP}
#' @param type a string specifying the type of plot to produce.
#' 'obs.dendrogram' produces the static cluster dendrogram for the observations;
#' 'var.dendrogram' produces the static cluster dendrogram for the variables;
#' 'heatmap' produces the static cluster heatmap with assiociated observation
#' and variable dendrograms;
#' and 'interactive' produces an interactive visualization of cluster
#' heatmap and its associated dendrograms.
#' @param ... Unused additional generic arguements
#' @param blwd a positive number. Line width on dendrograms.
#' @param lcex a positive number. Label size on dendrograms.
#' @param cexRow row label size
#' @param cexCol column label size
#' @export
#' @importFrom shiny shinyApp
#' @importFrom shiny fluidPage
#' @importFrom shiny titlePanel
#' @importFrom shiny tabsetPanel
#' @importFrom shiny fluidRow
#' @importFrom shiny animationOptions
#' @importFrom shiny tags
#' @importFrom shiny column
#' @importFrom shiny plotOutput
#' @importFrom shiny renderPlot
#' @importFrom shiny sliderInput
#' @importFrom stats as.dendrogram
#' @importFrom stats as.hclust
#' @importFrom stats quantile
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices adjustcolor
#' @examples
#' \dontrun{
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X[1:10,1:4]
#' cbass.fit <- CARP(
#'     X=Xdat,
#'     obs.labels=presidential_speech$labels[1:10])
#' plot(cbass.fit,type='interactive')
#' }
plot.CBASS <- function(
  x,
  type=c('obs.dendrogram','var.dendrogram','heatmap','interactive'),
  blwd=2,
  lcex=.6,
  cexRow=1,
  cexCol=1,
  ...){
  type = match.arg(type)
  switch(
    type,
    obs.dendrogram={
      x$cbass.dend.obs %>%
        stats::as.dendrogram() %>%
        dendextend::set("branches_lwd",blwd) %>%
        dendextend::set("labels_cex",lcex) %>%
        plot(ylab='Amount of Regularization')
    },
    var.dendrogram={
      x$cbass.dend.var %>%
        stats::as.dendrogram() %>%
        dendextend::set("branches_lwd",blwd) %>%
        dendextend::set("labels_cex",lcex) %>%
        plot(ylab='Amount of Regularization')
    },
    heatmap={
      if(x$X.center.global){
        X <- x$X
        X <- X - mean(X)
        X <- t(X)
      }else{
        X <- t(x$X.orig)
      }
      rownames(X) <- x$var.labels
      colnames(X) <- x$obs.labels
      nbreaks <- 50
      quant.probs <- seq(0,1,length.out = nbreaks)
      breaks <- unique(stats::quantile(X[TRUE],probs = quant.probs))
      nbreaks <- length(breaks)
      heatcols <- grDevices::colorRampPalette(c("blue","yellow"))(nbreaks - 1)

      my.cols <- grDevices::adjustcolor(c('black','grey'),alpha.f = .3)
      my.heatmap.2(x=X,
                   scale='none',
                   Colv=stats::as.dendrogram(x$cbass.dend.obs),
                   Rowv = stats::as.dendrogram(x$cbass.dend.var),
                   trace='none',
                   density.info = 'none',
                   key=FALSE,
                   breaks = breaks,
                   col=heatcols,
                   symkey = F,
                   Row.hclust = x$cbass.dend.var %>% stats::as.hclust(),
                   Col.hclust = x$cbass.dend.obs %>% stats::as.hclust(),
                   k.col=x$n.obs,
                   k.row=x$p.vars,
                   my.col.vec = my.cols,
                   cexRow = 1.5,
                   cexCol = 1.5,
                   margins = c(14,8))
    },
    interactive={
      shiny::shinyApp(
        ui=shiny::fluidPage(
          shiny::tags$style(type="text/css",
                     ".recalculating { opacity: 1.0; }"
          ),


          shiny::titlePanel("BiClustering"),
          shiny::tabsetPanel(
            shiny::tabPanel("Heatmap",
                     shiny::fluidRow(
                       shiny::column(3,
                              shiny::sliderInput(
                                "regcent",
                                "Amount of Regularization",
                                min = 0,
                                max = 1,
                                value = .5,
                                step=.03,animate = shiny::animationOptions(interval=300,loop=T)
                              )
                       ),
                       shiny::column(9,
                              shiny::plotOutput("heatmap",height = "900px",width = '1200px')
                       )
                     )
            )
          )
        ),
        server = function(input, output) {
          if(x$X.center.global){
            X.heat <- x$X
            X.heat <- X.heat - mean(X.heat)
            X.heat <- t(X.heat)
            X <- x$X
            X <- X - mean(X)
            X <- t(X)
          }else{
            X.heat <- t(x$X)
            X <- t(x$X)
          }
          colnames(X.heat) <- x$obs.labels
          rownames(X.heat) <- x$var.labels
          x$cbass.sol.path$lambda.path %>% as.vector() -> lam.seq
          lam.prop.seq <- lam.seq / max(lam.seq)
          nbreaks <- 50
          quant.probs <- seq(0,1,length.out = nbreaks)
          breaks <- unique(stats::quantile(X[TRUE],probs = quant.probs))
          nbreaks <- length(breaks)
          heatcols <- grDevices::colorRampPalette(c("blue","yellow"))(nbreaks - 1)
          my.cols <- grDevices::adjustcolor(c('black','grey'),alpha.f = .3)
          output$heatmap <- shiny::renderPlot({
            plt.iter <- which.min(abs(input$regcent - lam.prop.seq))
            # find lambda at iter
            cur.lam <- x$cbass.sol.path$lambda.path[plt.iter]
            cur.lam
            # find lambda closest in column path
            cur.col.lam.ind <- which.min(abs(x$cbass.cluster.path.obs$lambda.path.inter - cur.lam))
            # find clustering solution in column path
            cur.col.clust.assignment <- x$cbass.cluster.path.obs$clust.path[[cur.col.lam.ind]]$membership
            cur.col.clust.labels <- unique(cur.col.clust.assignment)
            cur.col.nclust <- length(cur.col.clust.labels)
            # find lambda closest in row path
            cur.row.lam.ind <- which.min(abs(x$cbass.cluster.path.var$lambda.path.inter - cur.lam))
            # find clustering solution in row path
            cur.row.clust.assignment <- x$cbass.cluster.path.var$clust.path[[cur.row.lam.ind]]$membership
            cur.row.clust.labels <- unique(cur.row.clust.assignment)
            cur.row.nclust <- length(cur.row.clust.labels)
            for(col.label.ind in seq_along(cur.col.clust.labels)){
              cur.col.label <- cur.col.clust.labels[col.label.ind]
              col.inds <- which(cur.col.clust.assignment == cur.col.label)
              for(row.label.ind in seq_along(cur.row.clust.labels)){
                cur.row.label <- cur.row.clust.labels[row.label.ind]
                row.inds <- which(cur.row.clust.assignment == cur.row.label)
                mean.value <- mean(X[row.inds,col.inds])
                X.heat[row.inds,col.inds] <- mean.value
              }
            }
            my.heatmap.2(x=X.heat,
                         scale='none',
                         Colv=stats::as.dendrogram(x$cbass.dend.obs),
                         Rowv = stats::as.dendrogram(x$cbass.dend.var),
                         trace='none',
                         density.info = 'none',
                         key=FALSE,
                         breaks = breaks,
                         col=heatcols,
                         symkey = F,
                         Row.hclust = x$cbass.dend.var %>% stats::as.hclust(),
                         Col.hclust = x$cbass.dend.obs %>% stats::as.hclust(),
                         k.col=cur.col.nclust,
                         k.row=cur.row.nclust,
                         my.col.vec = my.cols,
                         cexRow = 1,
                         cexCol = 1,
                         margins = c(10,10))

          })


    }
      )
    }
      )
}
