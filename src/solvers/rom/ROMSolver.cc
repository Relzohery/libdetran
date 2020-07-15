/*
 * ROM_Manager.cc
 *
 *  Created on: Jul 8, 2020
 *      Author: rabab
 */

#include "ROMSolver.hh"

template <class D>
ROMSolver<D>::ROMSolver(SP_input inp, SP_mesh mesh, SP_material mat, std::string operator_type)
:d_input(inp),
 d_mesh(mesh),
 d_mat(mat),
 d_operator_type(operator_type)

{
 Require(inp);
 Require(mat);
 Require(mesh);

 ROMSolver::Set_FullOperators();
}

template <class D>
void ROMSolver<D>::Set_FullOperators()
{
	if (d_operator_type == "diffusion")

	{
	  SP_lossoperator A (new DiffusionLossOperator(d_input, d_mat, d_mesh, false, 0.0, false, 1.0));
	  SP_gainoperator B (new DiffusionGainOperator(d_input, d_mat, d_mesh, false));

	  d_A = A;
	  d_B = B;
	}

	else if (d_operator_type == "EnergyDependent")
	{
	  d_input->put<std::string>("outer_solver", "GMRES");
	  d_input->put<std::string>("eigen_solver", "GD");
	  SP_mg_solver mg_solver_ED;
	  mg_solver_ED = new FixedSourceManager<_1D>(d_input, d_mat, d_mesh, false, true);
	  mg_solver_ED->setup();
	  mg_solver_ED->set_solver();
	  d_A = new LHS_Operator_T(mg_solver_ED);
	  //d_A->compute_explicit("EnergyDependent");
	}

	else if (d_operator_type == "EnergyIndependent")
	{
	 SP_mg_solver mg_solver;
	 mg_solver = new FixedSourceManager<_1D>(d_input, d_mat, d_mesh, false, true);
	 mg_solver->setup();
	 mg_solver->set_solver();
	 d_A = new Operator_T(mg_solver);
	 //d_A->compute_explicit("/home/rabab/Desktop/EnergyIndependent");
	}

}

template <class D>
void ROMSolver<D>::Solve(SP_matrix d_U, SP_vector sol)
{
	detran::OperatorProjection P(1);

	P.SetOperators(d_A, d_U);

    int d_r = d_U->number_columns();

	SP_matrix Ar;
	Ar = new callow::MatrixDense(d_r, d_r);
	P.Project(Ar);
	//Ar->print_matlab("Ar.txt");

	SP_eigensolver eigensolver;
	eigensolver = Creator_T::Create(d_input);

	if (d_operator_type == "diffusion")
	{
	  P.SetOperators(d_B, d_U);
      SP_matrix Br;
	  Br = new callow::MatrixDense(d_r, d_r);
	  P.Project(Br);
	  //Br->print_matlab("Br.txt");
	  eigensolver->set_operators(Br, Ar);
	}

	else eigensolver->set_operators(Ar);

    SP_vector x_rom;
    SP_vector x1;
    x_rom = new callow::Vector(d_r, 0.0);
    x1 = new callow::Vector(d_r, 1.0);
    eigensolver->solve(x_rom, x1);

    d_keff =  eigensolver->eigenvalue();

    // reconstruct
    int d_n = d_U->number_rows();

    x_fom = new callow::Vector(d_n, 0.0);
    d_U->multiply(x_rom, sol);
    }


//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

SOLVERS_INSTANTIATE_EXPORT(ROMSolver<_1D>)
SOLVERS_INSTANTIATE_EXPORT(ROMSolver<_2D>)
SOLVERS_INSTANTIATE_EXPORT(ROMSolver<_3D>)
