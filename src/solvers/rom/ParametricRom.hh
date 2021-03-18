/*
 * ParametricRom.h
 *
 *  Created on: Jan 3, 2021
 *      Author: rabab
 */

#ifndef SOURCE_SRC_SOLVERS_ROM_PARAMETRICROM_HH_
#define SOURCE_SRC_SOLVERS_ROM_PARAMETRICROM_HH_

#include "callow/solver/LinearSolverCreator.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/solver/EigenSolverCreator.hh"
#include "OperatorProjection.hh"
#include "offline_stage.hh"
#include "ROMSolver.cc"
#include "DEIM.hh"

using namespace callow;

namespace detran
{
template <class D>
class ParametricRom
{
public:
	typedef detran_utilities::InputDB::SP_input       SP_input;
	typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
	typedef detran_material::Material::SP_material    SP_material;
	typedef Matrix::SP_matrix                         SP_matrix;
	typedef MatrixDense::SP_matrix                    SP_matrixDense;
	typedef std::vector<SP_matrixDense>               vec_matrix;
	typedef Vector::SP_vector                         SP_vector;
	typedef ROMSolver<D>                              Base;
	typedef detran::DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
	typedef detran::DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
	typedef LinearSolverCreator::SP_solver            SP_solver;
	typedef callow::EigenSolver::SP_solver            SP_eigensolver;
	typedef callow::EigenSolverCreator                Creator_T;

	ParametricRom(SP_matrixDense Uf, SP_matrixDense UL, SP_matrixDense UR,
			      SP_input input, SP_mesh mesh, SP_material mat, int r_deim);

	void Set_FullOperators();

	void MDEIM_R();

	void MDEIM_L();
   // offline
    void reduced_matrices();

    void compute_coefficients();

    void construct_Operator();

    void solve();

private:

	SP_mesh d_mesh;
	SP_material d_mat;
	SP_input d_input;
	SP_matrixDense d_Uf;
	SP_matrixDense d_UL;
	SP_matrixDense d_UR;
	SP_matrixDense d_ULL;
	SP_matrixDense d_URR;
	int* d_Il;
	int* d_IR;
	int d_r;

    SP_matrix d_A;
	/// Operator "B" in "Ax = \lambda B"
	SP_matrix d_B;

	vec_matrix M_L;
	vec_matrix M_R;

	SP_vector a_L;
	SP_vector a_R;
	SP_solver d_solver;

	SP_matrixDense AR;
	SP_matrixDense AL;
};
};

#endif /* SOURCE_SRC_SOLVERS_ROM_PARAMETRICROM_HH_ */
