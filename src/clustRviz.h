#define ARMA_USE_SUPERLU 0

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#define ARMA_64BIT_WORD


// From https://stackoverflow.com/a/26268017
template <typename T_full, typename T_ind>
T_full extract(const T_full& full, const T_ind& ind){
  Eigen::Index num_indices = ind.innerSize();
  T_full target(num_indices);

  for(Eigen::Index i = 0; i < num_indices; i++){
    target[i] = full[ind[i]];
  }

  return target;
}
