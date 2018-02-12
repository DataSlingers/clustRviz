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
#' @export
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
#' @param carp.fit a CARP object returned by \code{CARP}
#' @export
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
  cat('Number of Observations:', carp.fit$n.obs,'\n')
  cat('Number of Variables:', carp.fit$p.vars,'\n')
  cat('Pre-processing:',preprocess.string[c(carp.fit$X.center,carp.fit$X.scale)],'\n')
  cat('Weights: RBF Kernel, phi =',carp.fit$phi, 'k =',carp.fit$k,'\n')
  cat('Algorithm:',alg.string,'\n')
  cat('Visualizations:',viz.string[c(carp.fit$static,carp.fit$static,carp.fit$interactive)],'\n')

  cat('Raw Data:\n')
  carp.fit$X[1:min(5,carp.fit$n.obs),1:min(5,carp.fit$p.vars)]

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
#' @param  carp.fit a CARP object returned by \code{CARP}
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
#' @param blwd a positive number. Line width on dendrograms.
#' @param lcex a positive number. Label size on dendrograms.
#' @export
#' @import shiny
#' @examples
#' library(clustRviz)
#' data("presidential_speech")
#' Xdat <- presidential_speech$X
#' carp.fit <- CARP(
#'     X=presidential_speech$X,
#'     obs.labels=presidential_speech$labels)
#' plot(carp.fit,type='interactive')
plot.CARP <- function(
  carp.fit,
  type='dendrogram',
  axis = c('PC1','PC2'),
  blwd=2,
  lcex=.6,
  percent=1,
  max.nclust=9,
  min.nclust=1){

  switch(
    type,
    dendrogram={
      carp.fit$carp.dend %>%
        as.dendrogram() %>%
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
      plot.frame <- carp.fit$carp.cluster.path.vis[,plot.cols]
      names(plot.frame)[1:2] <- c('V1','V2')
      plot.frame %>%
        filter(LambdaPercent <= percent) %>%
        filter(Iter > carp.fit$burn.in) %>%
        ggplot(aes(x=V1,y=V2,group=Obs)) +
        geom_path(
          aes(x=V1,y=V2),
          linejoin = 'round',
          color='red',
          size=1
        )  +
        geom_point(
          aes(x=V1,y=V2),
          data=plot.frame %>% filter(Iter==1),
          color='black',
          size=I(2)
        ) +
        ggrepel::geom_text_repel(
          aes(x=V1,y=V2,label=ObsLabel),
          size=I(3),
          data=plot.frame %>% filter(Iter == 1)
        ) +
        guides(color=FALSE,size=FALSE) +
        theme(axis.title = element_text(size=15)) +
        theme(axis.text = element_text(size=10)) +
        xlab(axis[1]) +
        ylab(axis[2])

    },
    interactive={
      shiny::shinyApp(
        ui=fluidPage(
          tags$style(type="text/css",
                     ".recalculating { opacity: 1.0; }"
          ),
          titlePanel("Clustering Example"),
          tabsetPanel(
            tabPanel("Movie",
                     fluidRow(
                       column(2,
                              sliderInput(
                                "regcent_movie",
                                "Amount of Regularization",
                                min = 0,
                                max = 1,
                                value = 0,
                                # step=.07,animate = animationOptions(interval=700,loop=T)
                                step=.03,animate = animationOptions(interval=1000,loop=T)
                              ),
                              uiOutput("choose_columns")
                       ),
                       column(5,
                              plotOutput("pcapathplot_movie",height = "700px")#,width = '900px')
                       ),
                       column(5,
                              plotOutput("dendplot_movie",height='700px')
                       )
                     )
            ),
            tabPanel("Static",
                     fluidRow(
                       column(2,
                              sliderInput("regcent_static",
                                          "Number of Clusters",
                                          min = min.nclust,
                                          max = max.nclust,
                                          value = max.nclust,
                                          step=1),
                              uiOutput("choose_columns_static")

                       ),
                       column(5,
                              plotOutput("pcapathplot_static",height = "700px")#,width = '900px')
                       ),
                       column(5,
                              plotOutput("dendplot_static",height='700px')
                       )
                     )
            )
          )
        ),
        server=function(input,output){

          # Drop-down selection box for which data set
          # Check boxes
          output$choose_columns <- renderUI({

            # Get the data set with the appropriate name
            colnames <- paste('PC',1:4,sep='')

            # Create the checkboxes and select them all by default
            checkboxGroupInput("columns", "Choose columns",
                               choices  = colnames,
                               selected = colnames[1:2])
          })
          # Check boxes
          output$choose_columns_static <- renderUI({
            # Get the data set with the appropriate name
            colnames <- paste('PC',1:4,sep='')

            # Create the checkboxes and select them all by default
            checkboxGroupInput("columns_static", "Choose columns",
                               choices  = colnames,
                               selected = colnames[1:2])
          })


          output$dendplot_movie <- renderPlot({
            carp.fit$carp.cluster.path.vis %>%
              filter(LambdaPercent <= input$regcent_movie)  %>%
              select(NCluster) %>%
              unlist() %>%
              unname() %>%
              min -> ncl
            carp.fit$carp.dend %>%
              as.dendrogram() %>%
              dendextend::set("branches_lwd",2) %>%
              dendextend::set("labels_cex",.6) %>%
              plot(ylab='Amount of Regularization',cex.lab=1.5)
            my.cols <- adjustcolor(c('grey','black'),alpha.f = .2)
            my.rect.hclust(carp.fit$carp.dend,k=ncl,border=2,my.col.vec=my.cols,lwd=3)


          })

          output$pcapathplot_movie <- renderPlot({
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
              carp.fit$carp.cluster.path.vis %>%
                filter(LambdaPercent <= input$regcent_movie) %>%
                select(Iter) %>%
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
            carp.fit$carp.cluster.path.vis %>%
              dplyr::select_(.dots=rename.list) -> carp.fit$carp.cluster.path.vis.rename

            carp.fit$carp.cluster.path.vis.rename %>%
              filter(Iter == cl.iter) %>%
              select(Obs,Cluster,Var1,Var2) %>%
              rename(
                MaxVar1 = Var1,
                MaxVar2 = Var2
              ) -> cl.assgn
            carp.fit$carp.cluster.path.vis.rename %>%
              filter(Iter <= cl.iter) %>%
              left_join(
                cl.assgn,
                by=c('Obs')
              ) -> tmp
            tmp %>%
              filter(Iter > min.iter ) %>%
              ggplot(aes(x=Var1,y=Var2,group=Obs)) +
              geom_path(
                aes(x=Var1,y=Var2),
                linejoin = 'round',
                color='red'
              ) +
              geom_point(
                aes(x=MaxVar1,y=MaxVar2),
                data = tmp %>% filter(Iter==1),
                color='red',
                size=I(4)
              ) +
              geom_point(
                aes(x=Var1,y=Var2),
                data=tmp %>% filter(Iter==1),
                color='black',
                size=I(4)
              ) +
              geom_text(
                aes(x=Var1,y=Var2,label=ObsLabel),
                size=I(6),
                data=tmp %>% filter(Iter == 1)
              ) +
              guides(color=FALSE,size=FALSE) +
              theme(axis.title = element_text(size=25)) +
              theme(axis.text = element_text(size=20)) +
              xlab(input$columns[1]) +
              ylab(input$columns[2])

          })
          output$dendplot_static <- renderPlot({
            carp.fit$carp.dend %>%
              as.dendrogram() %>%
              dendextend::set("branches_lwd",2) %>%
              dendextend::set("labels_cex",.6) %>%
              plot(ylab='Amount of Regularization',cex.lab=1.5)
            my.cols <- adjustcolor(RColorBrewer::brewer.pal(n=input$regcent_static,'Set1'),alpha.f=.2)
            my.rect.hclust(carp.fit$carp.dend,k=input$regcent_static,border=2,my.col.vec=my.cols,lwd=3)
          })
          output$pcapathplot_static <- renderPlot({
            ncl <- input$regcent_static
            my.cols <- adjustcolor(RColorBrewer::brewer.pal(n=ncl,'Set1'))[order(unique(cutree(carp.fit$carp.dend,k=ncl)[carp.fit$carp.dend$order]))]
            carp.fit$carp.cluster.path.vis %>%
              distinct(Iter,NCluster) %>%
              filter(NCluster == ncl) %>%
              select(Iter) %>%
              unlist() %>%
              unname() %>%
              median() -> cl.iter
            cl.iter <- floor(cl.iter)
            rename.list <- list(Obs = 'Obs',
                                Cluster = 'Cluster',
                                Iter = 'Iter',
                                ObsLabel = 'ObsLabel',
                                NCluster = 'NCluster',
                                Var1 = input$columns_static[1],
                                Var2 = input$columns_static[2])
            carp.fit$carp.cluster.path.vis %>%
              select_(.dots=rename.list) -> carp.cluster.path.vis.rename

            carp.cluster.path.vis.rename %>%
              filter(Iter == cl.iter) %>%
              select(Obs,Cluster,Var1,Var2) %>%
              rename(
                MaxVar1 = Var1,
                MaxVar2 = Var2,
                PlotCluster=Cluster
              ) %>%
              mutate(
                PlotCluster = as.factor(PlotCluster)
              ) -> cl.assgn

            carp.cluster.path.vis.rename %>%
              filter(Iter <= cl.iter) %>%
              left_join(
                cl.assgn,
                by=c('Obs')
              ) -> tmp
            tmp %>%
              filter(Iter > 50 ) %>%
              ggplot(aes(x=Var1,y=Var2,group=Obs)) +
              geom_path(
                aes(x=Var1,y=Var2,color=PlotCluster),
                linejoin = 'round'
              ) +
              geom_point(
                aes(x=MaxVar1,y=MaxVar2,color=PlotCluster),
                data = tmp %>% filter(Iter==1),
                size=I(4)
              ) +
              geom_point(
                aes(x=Var1,y=Var2),
                data=tmp %>% filter(Iter==1),
                color='black',
                size=I(4)
              ) +
              geom_text(
                aes(x=Var1,y=Var2,label=ObsLabel),
                size=I(6),
                data=tmp %>% filter(Iter == 1)
              ) +
              scale_color_manual(values=my.cols)+
              guides(color=FALSE,size=FALSE) +
              theme(axis.title = element_text(size=25)) +
              theme(axis.text = element_text(size=20)) +
              xlab(input$columns_static[1]) +
              ylab(input$columns_static[2])
          })



        }
        # End Shiny App
      )

    }
  )

}
