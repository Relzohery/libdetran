/*
 * DEIM.hh
 *
 *  Created on: Jan 1, 2021
 *      Author: rabab
 */

#ifndef SOURCE_SRC_SOLVERS_ROM_DEIM_HH_
#define SOURCE_SRC_SOLVERS_ROM_DEIM_HH_

#include "solvers/solvers_export.hh"
#include "callow/solver/LinearSolverCreator.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/matrix/MatrixShell.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/vector/Vector.hh"
#include "OperatorProjection.hh"
#include "utilities/InputDB.hh"
#include <iostream>

using namespace callow;

namespace detran
{
class DEIM
{
  public:
    typedef MatrixDense::SP_matrix  SP_matrixDense;
	typedef MatrixBase::SP_matrix   SP_matrix;
	typedef LinearSolverCreator::SP_solver            SP_solver;
	typedef callow::Vector::SP_vector                 SP_vector;
	typedef detran_utilities::InputDB::SP_input       SP_input;

	DEIM(SP_matrixDense U, int r);

	virtual ~DEIM(){}

	void Search();

	/// getter of the interpolation indices
	int* interpolation_indices()
	{
      return d_l;
	}

	SP_matrixDense ReducedBasis()
	{
	  return Ur;
	}

  private:
    /// DEIM rank
    unsigned int d_r;
    /// problem size
    unsigned int d_n;
    /// interpolation indices
    int* d_l;
    /// operator basis
    SP_matrixDense d_U;
    /// reduced basis
    SP_matrixDense d_Um;
    /// linear solver
    SP_solver d_solver;
    ///
    SP_matrixDense Ur;

    LinearSolverCreator::SP_db db;
    LinearSolver::SP_db get_db()
    {
      LinearSolver::SP_db p(new detran_utilities::InputDB("callow_db"));

	  p->put<double>("linear_solver_atol",                 1e-16);
	  p->put<double>("linear_solver_rtol",                 1e-15);
	  p->put<std::string>("linear_solver_type", "gmres");

	  //p->put<std::string>("linear_solver_type", "petsc");
	  //p->put<std::string>("pc_type", "petsc_pc");
	  //p->put<std::string>("petsc_pc_type", "lu");
      p->put<int>("linear_solver_maxit", 50);
	  p->put<int>("linear_solver_gmres_restart", 16);
	  p->put<int>("linear_solver_maxit",                   2000);
	  p->put<int>("linear_solver_gmres_restart",           30);
	  p->put<int>("linear_solver_monitor_level",           0);

	  return p;
 }
    void findIndex(int col);
};
} /// end of namespace

#endif /* SOURCE_SRC_SOLVERS_ROM_DEIM_HH_ */
