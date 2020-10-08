#ifndef CLUSTRVIZ_TROUT_H
#define CLUSTRVIZ_TROUT_H 1

#include "clustRviz_base.h"
#include "clustRviz_logging.h"
#include "status.h"

template <class DataType>
class UnivariateTroutClusteringSkeleton : public ConvexClusteringSkeleton<DataType> {
  using MatrixXt = Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorXt = Eigen::Matrix<DataType, Eigen::Dynamic, 1>;

  // We need to inherit the constructor explicitly for reasons I don't entirely understand
  // or else this fails with a confusing error message which essentially says the
  // constructor call (with 7 arguments) doesn't match the default constructor
  using ConvexClusteringSkeleton<DataType>::ConvexClusteringSkeleton;

  // TODO - Allow flexible specification of phase and amplitude alignment
  // TODO - Allow U-step convergence parameters to be tweaked
  inline void u_step(){
    // For now, we do a whole loop, but at some point might relax this to a single step if theory supports...

    // Note that "this->" is required for member access since we're inheriting from a
    // templated base class -- https://stackoverflow.com/a/1121016
    Eigen::Index k_inner = 0;
    MatrixXt U_old_inner = this->U;
    const MatrixXt A = this->rho * this->D.transpose() * (this->V - this->Z);

    do{
      k_inner += 1;
      U_old_inner = this->U;
      MatrixXt X_aligned = align_phase(this->X, this->U);

      this->U = this->u_step_solver.solve(X_aligned + A);

    } while ((scaled_squared_norm(this->U - U_old_inner) > 1e-6) & (k_inner < 1000));

    ClustRVizLogger::debug("U = ") << this->U;
  }
};

#endif
