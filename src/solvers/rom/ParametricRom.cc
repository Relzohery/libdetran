/*
 * ParametricRom.cpp
 *
 *  Created on: Jan 3, 2021
 *      Author: rabab
 */

#include "ParametricRom.hh"
namespace detran
{
template <class D>
ParametricRom<D>::ParametricRom(SP_matrixDense Uf, SP_matrixDense UL, SP_matrixDense UR,
	      SP_input input, SP_mesh mesh, SP_material mat, int r_deim)
:d_input(input),
 d_mesh(mesh),
 d_mat(mat),
 d_UL(UL),
 d_UR(UR),
 d_Uf(Uf),
 d_r(r_deim)
{
 Require(input);
 Require(mat);
 Require(mesh);

 MDEIM_L();

 MDEIM_R();

 Set_FullOperators();

 reduced_matrices();

 compute_coefficients();

 construct_Operator();

 //solve();
}

template <class D>
void ParametricRom<D>::Set_FullOperators()
{
  SP_lossoperator A (new DiffusionLossOperator(d_input, d_mat, d_mesh, false, 0.0, false, 1.0, false));

  SP_gainoperator B (new DiffusionGainOperator(d_input, d_mat, d_mesh, false));

   d_A = A;
   d_A->number_columns();
   d_B = B;
}
template <class D>
void ParametricRom<D>::MDEIM_L()
{
  // the left operator
  DEIM DEIM_L(d_UL, d_r);
  DEIM_L.Search();
  d_Il = DEIM_L.interpolation_indices();
  d_ULL = new callow::MatrixDense(d_r, d_r);
  d_ULL = DEIM_L.ReducedBasis();
}

template <class D>
void ParametricRom<D>::MDEIM_R()
{
  // the right operator
  DEIM DEIM_R(d_UR, d_r);
  DEIM_R.Search();
  d_IR = DEIM_R.interpolation_indices();
  d_URR = new callow::MatrixDense(d_r, d_r);
  d_URR = DEIM_R.ReducedBasis();

}

template <class D>
void ParametricRom<D>::reduced_matrices()
{
  offline_stage off_L(d_A, d_Uf, d_UL, d_r);
  M_L = off_L.Decompositon();
  //
  offline_stage off_R(d_B, d_Uf, d_UR, d_r);
  M_R = off_R.Decompositon();

}

///////////// online /////////////////

template <class D>
void ParametricRom<D>::compute_coefficients()
{
  LinearSolverCreator::SP_db db;
  LinearSolver::SP_db p(new detran_utilities::InputDB("callow_db"));
  p->put<int>("linear_solver_maxit",   1000);
  p->put<double>("linear_solver_atol", 1e-13);
  p->put<double>("linear_solver_rtol", 1e-13);
  p->put<int>("linear_solver_monitor_level", 0);
  p->put<int>("linear_solver_monitor_diverge", 0);
  p->put<std::string>("linear_solver_type", "petsc");
  p->put<std::string>("pc_type", "lu");
  d_solver = LinearSolverCreator::Create(db);

  Vector b(d_r, 0.0);
  a_L = new Vector (d_r, 0.0);
  double* nnz_values = d_A->values();
  for (int i=0; i<d_r; i++)
  {
    b[i] = nnz_values[d_Il[i]];
    std::cout << b[i] << "   \n";
  }

  b.print_matlab("b.txt");
  d_ULL->print_matlab("A_DEIM.txt");

  d_solver->set_operators(d_ULL);
  d_solver->solve(b, *a_L);

  a_R = new Vector (d_r, 0.0);
  double* nnz_values_R = d_B->values();
  for (int i=0; i<d_r; i++)
  {
    b[i] = nnz_values_R[d_IR[i]];
  }

  b.print_matlab("b_G.txt");
  d_ULL->print_matlab("A_DEIM_G.txt");
  d_solver->set_operators(d_URR);
  d_solver->solve(b, *a_R);
}

//
template <class D>
void ParametricRom<D>::construct_Operator()
{
  double v;
//  // construct the Left operator
  AL = new MatrixDense(14, 14);
  for (int r=0; r< d_r; r++)
  {
   v = 0;
   for (int i=0; i<14; i++)
   {
    for (int j=0; j<14; j++)
    {
     v = (*a_L)[r]*(*M_L[r])(i, j);
     AL->insert(i, j, v, 1);
    }
   }
 }

  AL->print_matlab("A_L.txt");

  AR = new MatrixDense(14, 14);
  for (int r=0; r< d_r; r++)
  {
    v = 0;
    for (int i=0; i<14; i++)
    {
      for (int j=0; j<14; j++)
    {
      v = (*a_R)[r]*(*M_R[r])(i, j);
      AR->insert(i, j, v, 1);
    }
   }
  }

  AR->print_matlab("A_R.txt");
}
//

template <class D>
void ParametricRom<D>::solve()
{

  AR->print_matlab("A_R.txt");
  AL->print_matlab("A_L.txt");
  SP_eigensolver eigensolver;
  eigensolver = Creator_T::Create(d_input);
  eigensolver->set_operators(AR, AL);
  SP_vector x_rom;
  SP_vector x;
  x_rom = new Vector(14, 0.0);
  x = new Vector(14, 1.0);
  eigensolver->solve(x_rom, x);
  double _keff =  eigensolver->eigenvalue();

  std::cout << M_L[0]->number_rows() << "\n";
  std::cout << d_r << "\n";
  std::cout << _keff << " &&&&&&&&&&&&\n";
  // reconstruct
  int d_n = d_Uf->number_rows();
  SP_vector sol;
  sol = new callow::Vector(d_n, 0.0);

  //d_Uf->multiply(*x_rom, *sol);
//  // correct the direction of the vector if negative
//  for (int i=0; i<d_n; i++)
//  {
//   if (((*sol)[i]) < 0)
//   {
//     (*sol)[i] *= -1;
//   }
//   }
//
}
//

template class ParametricRom<_1D>;
template class ParametricRom<_2D>;
template class ParametricRom<_3D>;
}
