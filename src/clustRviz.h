#include "clustRviz_base.h"
#include "clustRviz_logging.h"
#include "clustering_impl.h"
#include "biclustering_impl.h"
#include "alg_reg_policies.h"
#include "optim_policies.h"

typedef AlgorithmicRegularizationFixedStepSizePolicy<ConvexClustering> CARP;
typedef AlgorithmicRegularizationBacktrackingPolicy<ConvexClustering> CARP_VIZ;
typedef AlgorithmicRegularizationFixedStepSizePolicy<ConvexBiClustering> CBASS;
typedef AlgorithmicRegularizationBacktrackingPolicy<ConvexBiClustering> CBASS_VIZ;
typedef ADMMPolicy<ConvexClustering> ConvexClusteringADMM;
typedef ADMMPolicy<ConvexBiClustering> ConvexBiClusteringADMM;
typedef BackTrackingADMMPolicy<ConvexClustering> ConvexClusteringADMM_VIZ;
typedef BackTrackingADMMPolicy<ConvexBiClustering> ConvexBiClusteringADMM_VIZ;
