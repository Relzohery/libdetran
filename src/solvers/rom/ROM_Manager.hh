/*
 * ROM_Solver.hh
 *
 *  Created on: Jul 8, 2020
 *      Author: rabab
 */

#ifndef SOLVERS_ROM_ROM_MANAGER_HH_
#define SOLVERS_ROM_ROM_MANAGER_HH_


#include "callow/vector/Vector.hh"
#include "callow/matrix/MatrixBase.hh"
#include "utilities/InputDB.hh"
#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "callow/solver/EigenSolverCreator.hh"
#include "FixedSourceManager.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/eigen/EnergyDependentEigenLHS.hh"
#include "solvers/eigen/EnergyIndependentEigenOperator.hh"
#include "ProjectedOperator.hh"

#include <string>

using namespace detran;

template <class D>
class ROM_Manager
{
public:
   typedef callow::EigenSolver::SP_solver            SP_eigensolver;
   typedef callow::EigenSolverCreator                Creator_T;
   typedef callow::MatrixBase::SP_matrix             SP_matrix;
   typedef callow::Vector::SP_vector                 SP_vector;
   typedef detran_utilities::InputDB::SP_input       SP_input;
   typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
   typedef detran_material::Material::SP_material        SP_material;
   typedef detran::DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
   typedef detran::DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
   typedef typename Eigensolver<D>::Fixed_T              Fixed_T;
   typedef typename Eigensolver<D>::SP_solver            SP_solver;
   typedef typename Fixed_T::SP_manager                  SP_mg_solver;
   typedef EnergyDependentEigenLHS<D>                LHS_Operator_T;
   typedef EnergyIndependentEigenOperator<D>             Operator_T;




   ROM_Manager(SP_input inp, SP_mesh Mesh, SP_material mat, std::string operator_type);

   void SetBasis();

   void Set_FullOperators();

   void SetOperators(SP_matrix A, SP_matrix B);

   void Solve(SP_matrix d_U);

   void ComputeROM_Error(std::string type);

   double keff()
   {
     return d_keff;
   }

  callow::Vector fom_state()
  {
	return x_fom;
  }

  callow::Vector rom_state()
  {
    return x_rom;
  }

private:
	SP_input d_input ;
	SP_mesh d_mesh;
	SP_material d_mat;
	SP_matrix  d_A ;
	SP_matrix d_B;

	std::string d_operator_type;
	double d_keff;
	callow::Vector x_rom;
	callow::Vector x_fom;
};

#endif /* SOLVERS_ROM_ROM_MANAGER_HH_ */
