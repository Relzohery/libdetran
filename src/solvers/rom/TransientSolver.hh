//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TransientSolver.hh
 *  @brief TransientSolver class definition.
 *  @note  Copyright(C) 2020 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef SOLVERS_ROM_TRANSIENTSOLVER_HH_
#define SOLVERS_ROM_TRANSIENTSOLVER_HH_

#include "callow/utils/Initialization.hh"
#include "callow/solver/LinearSolverCreator.hh"
#include "callow/solver/LinearSolver.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "kinetics/TimeDependentMaterial.hh"
#include "kinetics/KineticsMaterial.hh"
#include "kinetics/MultiPhysics.hh"
#include "time/TimeStepper.hh"
#include "callow/matrix/MatrixShell.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/vector/Vector.hh"
#include "transport/State.hh"
#include "kinetics/Precursors.hh"
#include "kinetics/SyntheticSource.hh"
#include "utilities/Definitions.hh"
#include "OperatorProjection.hh"
#include "transport/State.hh"
#include "KineticMatrices.hh"
#include "offline_stage.hh"
#include "DEIM.hh"

using namespace detran;
using namespace callow;


class TransientSolver
{

public:
  typedef callow::MatrixDense::SP_matrix            SP_matrix;
  typedef callow::Vector::SP_vector                 SP_vector;
  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef detran_utilities::vec_int                 vec_int;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef TimeDependentMaterial::SP_material        SP_material;
  typedef Precursors::SP_precursors                 SP_precursors;
  typedef FissionSource::SP_fissionsource           SP_fissionsource;
  typedef State::SP_state                           SP_state;
  typedef DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
  typedef DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
  typedef LinearSolverCreator::SP_solver            SP_solver;
  typedef callow::Matrix::SP_matrix                 SP_matrix_Base;
  typedef std::vector<SP_matrix>                    vec_matrix;
  typedef detran_utilities::SP<MultiPhysics>        SP_multiphysics;
  typedef std::vector<callow::Vector>               vec_flux;

  typedef void (*multiphysics_pointer)
                (void*, SP_vector , double, double, vec_matrix, SP_matrix);
  typedef std::vector<SP_multiphysics>              vec_multiphysics;


  TransientSolver(SP_input inp, SP_mesh mesh, SP_material material, SP_matrix flux_basis, SP_matrix precursors_basis);

  void Solve(SP_state initial_state);

  /// set DEIM basis
  void set_DEIM(SP_matrix U_deim);

  void set_multiphysics(SP_multiphysics ic,
                      multiphysics_pointer update_multiphysics_rhs,
					  SP_matrix Temp_State_basis,
                      void* multiphysics_data = NULL);
  /// flux getter
  SP_matrix flux() {return d_flux;};
  /// precursors getter
  SP_matrix precursors() {return d_precursors;};
  /// power getter
  SP_vector power() {return d_power;};


private:
  /// State vector
  SP_state d_state;
  /// Fission source
  SP_fissionsource d_fissionsource;
  /// Precursors production matrix
  SP_matrix d_precursors_production;
  /// Precursors decsy matrix
  SP_matrix d_precursors_decay;
  /// Delayed neutron production
  SP_matrix d_delayed_production;
  /// Gain matrix
  SP_matrix_Base d_G;
  /// Loss matrix
  SP_matrix_Base d_L;
  /// operator
  SP_matrix d_A;
  /// implicit Euler operator
  SP_matrix d_A_;
  SP_matrix d_LF;
  SP_matrix d_Lr;
  SP_matrix d_Gr;
  /// Reduced solution
  SP_matrix d_sols;
  /// Flux solution matrix at all time steps
  SP_matrix d_flux;
  /// Precursors solution matrix at all time steps
  SP_matrix d_precursors;
  /// Reduced flux vector
  SP_matrix d_flux_r;
  /// Reduced precursors vector
  SP_matrix d_precursors_r;
  double* vectorized_matrix;

  /// Time step size
  double d_dt;
  /// Final time
  double  d_final_time;
  /// Number of time steps
  double d_number_steps;

  /// Mesh
  SP_mesh d_mesh;
  /// Input
  SP_input d_inp;
  /// Solver setting
  SP_input db;
  /// Material
  SP_material d_material;
  /// Flux basis
  SP_matrix d_flux_basis;
  /// Precursors basis
  SP_matrix d_precursors_basis;
  /// DEIM basis
  SP_matrix d_deim_basis;
  /// Precursor vector in previous time step
  SP_vector d_P0;
  /// Flux vector
  SP_vector d_phi;
  /// Current vector of the reduced flux and precursors
  SP_vector d_sol_r;
  /// Vector of the reduced flux and precursors at the previous time step
  SP_vector d_sol0_r;
  /// vector of fixed point iteration
  SP_vector d_x0;
  /// The projected vector of the precursors
  SP_vector d_P_r;
  /// The projected vector of the initial precursors
  SP_vector d_P0_r;
  /// The projected vector of the initial flux
  SP_vector d_phi0_r;
  /// The projected vector of the flux
  SP_vector d_phi_r;
  /// reduced DEIM mats
  vec_matrix M_L;
  /// reduced deim basis
  SP_matrix Ur_deim;
  /// coefficients of the decomposed matrix
  SP_vector d_x_deim;
  /// selected elements of the vectorized matrix
  SP_vector d_b_deim;
  /// power
  SP_vector d_power;
  /// DEIM interpolation indices
  int* l;
  ///
  SP_matrix d_L_deim;
  ///
  bool deim_flag = false;
  /// Number of cells
  int d_num_cells;
  /// number of energy groups
  int d_number_groups;
  /// Number of precursors group
  int d_precursors_group ;
  /// Flux rank
  int d_rf;
  /// Precursors rank
  int d_rc;

  int d_n;
  /// deim rank
  int r_deim;
  /// Linear solver
  SP_solver d_solver;
  /// linear solver DEIM
  SP_solver d_solver_deim;

  double d_tolerance;
  /// Maximum nonlinear iterations
  size_t d_maximum_iterations;
  /// Compute the initial precursors concentration
  void initialize_precursors();
  /// Project the initial flux and precursors on space of the reduced basis
  void ProjectInitial();
  /// Construct matrix
  void Construct_Operator(double t, double dt);
  /// Update the operator 
  void Refersh_Operator();
  /// Reconstruct the full order solution
  void reconstruct(int i);
  /// Working multiphysics vector
  SP_multiphysics d_multiphysics;
  /// Multiphysics callback
  multiphysics_pointer d_update_multiphysics_rhs;
  /// Multiphysics data
  void* d_multiphysics_data;
  /// Vector of previous physics iterates
  vec_multiphysics d_vec_multiphysics;
  /// Previous multiphysics iterate
  SP_multiphysics d_multiphysics_0;

  vec_matrix Temp_State_basis;

  SP_matrix multiphysics_basis;

  SP_matrix XS;
  void step(int, double);

  /// DEIM offline stage
  void DEIM_offline();

  /// DEIM online stage
  void DEIM_online();
  /// update multiphysics at time step
  void update_multiphysics(const double t, const double dt, const size_t order, vec_matrix, SP_matrix);
  /// check flux convergence at each iteration
  bool check_convergence();


};

#endif /* SOLVERS_ROM_TRANSIENTSOLVER_HH_ */
