#include "clustRviz_base.h"
#include "clustRviz_logging.h"
#include "clustRviz_utils.h"
#include "clustering_impl.h"
#include "biclustering_impl.h"
#include "trout_impl.h"
#include "alg_reg_policies.h"
#include "optim_policies.h"

typedef ConvexClusteringSkeleton<double> ConvexClustering;
typedef UnivariateTroutClusteringSkeleton<std::complex<double> > TroutClustering;

typedef AlgorithmicRegularizationFixedStepSizePolicy<ConvexClustering> CARP;
typedef AlgorithmicRegularizationBacktrackingPolicy<ConvexClustering> CARP_VIZ;
typedef AlgorithmicRegularizationFixedStepSizePolicy<ConvexBiClustering> CBASS;
typedef AlgorithmicRegularizationBacktrackingPolicy<ConvexBiClustering> CBASS_VIZ;
typedef AlgorithmicRegularizationFixedStepSizePolicy<TroutClustering> TROUT;
typedef AlgorithmicRegularizationBacktrackingPolicy<TroutClustering> TROUT_VIZ;

typedef ADMMPolicy<ConvexClustering> ConvexClusteringADMM;
typedef ADMMPolicy<ConvexBiClustering> ConvexBiClusteringADMM;
typedef ADMMPolicy<TroutClustering> TroutClusteringADMM;
typedef BackTrackingADMMPolicy<ConvexClustering> ConvexClusteringADMM_VIZ;
typedef BackTrackingADMMPolicy<ConvexBiClustering> ConvexBiClusteringADMM_VIZ;
typedef BackTrackingADMMPolicy<TroutClustering> TroutClusteringADMM_VIZ;
typedef UserGridADMMPolicy<ConvexClustering> UserGridConvexClusteringADMM;
typedef UserGridADMMPolicy<ConvexBiClustering> UserGridConvexBiClusteringADMM;
typedef UserGridADMMPolicy<TroutClustering> UserGridTroutClusteringADMM;
