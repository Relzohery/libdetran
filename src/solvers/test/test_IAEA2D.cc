/*
 * test_IAEA2D.cc
 *
 *  Created on: Mar 22, 2021
 *      Author: rabab
 */


#define TEST_LIST \
        FUNC(test_IAEA2D)\
		FUNC(test_IAEA2D_ROM)\
		FUNC(test_DEIM_snapshots)\
		FUNC(test_IAEA_DEIM)


#include "TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "utilities/MathUtilities.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "solvers/EigenvalueManager.hh"
#include "callow/vector/Vector.hh"
#include "solvers/rom/ROMSolver.hh"
#include "solvers/rom/ROMBasis.hh"
#include "Mesh2D.hh"
#include "DEIM.hh"
#include "offline_stage.hh"
#include <iostream>
#include <fstream>

using namespace detran;
using namespace detran_test;
using namespace detran_utilities;
using namespace detran_material;
using namespace detran_geometry;
using namespace callow;
using namespace std;
using std::cout;
using std::endl;

typedef Mesh::SP_mesh SP_mesh;


typedef callow::MatrixDense::SP_matrix SP_matrix;
typedef callow::Vector::SP_vector      SP_vector;



int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}


double** get_samples(unsigned m, unsigned n)
{
  double** data = 0;
  data =  new double* [m];

  std::ifstream input("../../../source/src/solvers/test/rom_basis/IAEA2D_samples_test.txt");
  for (int i = 0; i < m; i++)
  {
	data[i] = new double[n];
    for (int j=0; j <n; j++)
	{
	  input >> data[i][j];
	  std::cout<< data[i][j]<<std::endl;
	}
  }
 return data;
}



InputDB::SP_input get_input()
 {
   InputDB::SP_input inp(new InputDB("IAEA-2D"));
   inp->put<std::string>("problem_type", "eigenvalue");
   inp->put<int>("number_groups",                      2);
   inp->put<int>("dimension",                          2);
   inp->put<std::string>("equation",                           "diffusion");
   inp->put<std::string>("bc_west",                            "reflect");
   inp->put<std::string>("bc_east",                            "vacuum");
   inp->put<std::string>("bc_south",                           "reflect");
   inp->put<std::string>("bc_north",                           "vacuum");
   inp->put<std::string>("inner_solver",           "SI");
   inp->put<int>("inner_max_iters",            1000);
   inp->put<int>("inner_max_iters",            1000);
   inp->put<double>("inner_tolerance",            1e-7);
   inp->put<int>("inner_print_level",          0);
   inp->put<std::string>("outer_solver",              "GS");
   inp->put<int>("outer_max_iters",            1000);
   inp->put<double>("outer_tolerance",            1e-7);
   inp->put<int>("inner_print_level",          0);
   inp->put<int>("outer_print_level",          0);
   //inp->put<std::string>("eigen_solver",       "arnoldi");
   inp->put<int>("eigen_max_iters",            200);
   inp->put<double>("eigen_solver_tol",            1e-14);
   inp->put<int>("quad_number_polar_octant",   16);

    InputDB::SP_input db(new InputDB("callow dp"));
    db->put<double>("linear_solver_atol",              1e-12);
    db->put<double>("linear_solver_rtol",              1e-12);
    db->put<std::string>("linear_solver_type",              "petsc");
    db->put<int>("linear_solver_maxit",             5000);
    db->put<int>("linear_solver_gmres_restart",     30);
    db->put<int>("linear_solver_monitor_level",     0);
    db->put<std::string>("pc_type",                         "petsc_pc");
    db->put<std::string>("petsc_pc_type",                   "lu");
    db->put<std::string>("eigen_solver_type",               "slepc");
    db->put<double>("eigen_solver_tol", 1e-14);
    db->put<int>("eigen_solver_monitor_level",      2);
    inp->put<InputDB::SP_input>("inner_solver_db",               db);
    inp->put<InputDB::SP_input>("inner_pc_db",                   db);
    inp->put<InputDB::SP_input>("outer_solver_db",               db);
    inp->put<InputDB::SP_input>("eigen_solver_db",               db);
    inp->put<InputDB::SP_input>("rom_solver_db",                 db);

   return inp;

 }

Material::SP_material get_mat(double* xs)
{
  Material::SP_material mat = Material::Create(4, 2, "IAEA_2D");
  // Material 0

  mat->set_sigma_t(0, 0, xs[0]);
  mat->set_sigma_t(0, 1, xs[1]);
  mat->set_sigma_s(0, 1, 0, xs[3]);
  mat->set_sigma_f(0, 1, xs[9]);
  mat->set_chi(0, 0, 1.0);
  mat->set_diff_coef(0, 0, xs[5]);
  mat->set_diff_coef(0, 1, xs[6]);
  mat->set_chi(0, 0, 1.0);
  mat->compute_sigma_a();

  //mat.set_diff_coef(0, vec_dbl([1.5, 0.4]))

  //Material 1
  mat->set_sigma_t(1, 0, xs[0]);
  mat->set_sigma_t(1, 1, xs[12]);
  mat->set_sigma_s(1, 1, 0, xs[3]);
  mat->set_sigma_f(1, 1, xs[9]);
  mat->set_diff_coef(1, 0, xs[5]);
  mat->set_diff_coef(1, 1, xs[6]);
  mat->set_chi(1, 0, 1.0);
  mat->compute_sigma_a();

  // Material 2
  mat->set_sigma_t(2, 0, xs[0]);
  mat->set_sigma_t(2, 1, xs[2]);
  mat->set_sigma_s(2, 1, 0, xs[3]);
  mat->set_sigma_f(2, 1,  xs[9]);
  mat->set_diff_coef(2, 0, xs[5]);
  mat->set_diff_coef(2, 1, xs[6]);
  mat->set_chi(2, 0, 1.0);
  mat->compute_sigma_a();

  // Material 3
  mat->set_sigma_t(3, 0, xs[10]);
  mat->set_sigma_t(3, 1, xs[11]);
  mat->set_sigma_s(3, 1, 0, xs[4]);
  mat->set_diff_coef(3, 0, xs[7] );
  mat->set_diff_coef(3, 1, xs[8]);
  mat->set_chi(3, 0, 1.0);
  mat->compute_sigma_a();
  mat->finalize();
  //mat->display();

 return mat;
}

Mesh2D::SP_mesh get_mesh(int fmm = 1, std::string id="core")
{
  Mesh2D::vec_dbl cm(10);
  cm[0] = 0.0;
  cm[1] = 10.0;
  cm[2] = 30.0;
  cm[3] = 50.0;
  cm[4] = 70.0;
  cm[5] = 90.0;
  cm[6] = 110.0;
  cm[7] = 130.0;
  cm[8] = 150.0;
  cm[9] = 170.0;

  Mesh2D::vec_int fm(9);
  for (int i=0; i < 9; i++) fm[i] = 10;

  int  mt_map[] = {2, 1, 1, 1, 2, 1, 1, 0, 3,
		           1, 1, 1, 1, 1, 1, 1, 0, 3,
		           1, 1, 1, 1, 1, 1, 0, 0, 3,
		           1, 1, 1, 1, 1, 1, 0, 3, 3,
		           2, 1, 1, 1, 2, 0, 0, 3, 3,
		           1, 1, 1, 1, 0, 0, 3, 3, 3,
		           1, 1, 0, 0, 0, 3, 3, 3, 3,
		           0, 0, 0, 3, 3, 3, 3, 3, 3,
		           3, 3, 3, 3, 3, 3, 3, 3, 3};

  Mesh2D::vec_int mt(9*9);
    for (int i = 0; i < 81; ++i) mt[i] = mt_map[i];

  //vec_int mat_map(9 * 9);
//  for (int i = 0; i < mt.size(); ++i)
//  mat_map[i] = mt[i];
  Mesh2D::SP_mesh mesh = Mesh2D::Create(fm, fm, cm, cm, mt);


return mesh;
}

//--------------------------------------------------------------------------//

int test_IAEA2D(int argc, char *argv[])
{
  time_t begin, end;
  time(&begin);

  InputDB::SP_input inp = get_input();

  //-------------------------------------------------------------------------//
  // MATERIAL
  //-------------------------------------------------------------------------//
  inp->get<std::string>("equation") = "diffusion";

  Mesh2D::SP_mesh mesh = get_mesh();

  MatrixDense phi_0(mesh->number_cells(), 100);
  MatrixDense phi_1(mesh->number_cells(), 100);
  Vector k(100);

  double** data = get_samples(100, 13);
  for (int i=0; i<100; i++)
  {
	Material::SP_material mat = get_mat(data[i]);
    EigenvalueManager<_2D> manager(inp, mat, mesh);
    manager.solve();
    State::SP_state ic = manager.state();

    std::cout << manager.state()->phi(0).size() << "\n";

    double* a ;
	double* b;
    a = new double [manager.state()->phi(0).size()];
    b = new double [manager.state()->phi(1).size()];

    for (int i=0; i< 8100; i++)
    {
      a[i] = manager.state()->phi(0)[i];
      b[i] = manager.state()->phi(1)[i];
    }
    phi_0.insert_col(i, a);
    phi_1.insert_col(i, b);
    k[i] = manager.state()->eigenvalue();

  }

  phi_0.print_matlab("IAEA_FOM_flux0.txt");
  phi_1.print_matlab("IAEA_FOM_flux1.txt");
  k.print_matlab("IAEA_FOM_eigenvalue.txt");

  return 0;
}


int test_IAEA2D_ROM(int argc, char *argv[])
{
  time_t begin, end;
  time(&begin);
  InputDB::SP_input inp = get_input();

  //-------------------------------------------------------------------------//
  // MATERIAL
  //-------------------------------------------------------------------------//
  inp->get<std::string>("equation") = "diffusion";
  inp->put<std::string>("operator", "diffusion");

  Mesh2D::SP_mesh mesh = get_mesh();

  int r = 50;
  int n = mesh->number_cells();
  // get the basis
  SP_matrix U;
  U = new callow::MatrixDense(2*n, 2*r);

  ROMBasis::GetBasis("../../../source/src/solvers/test/rom_basis/IAEA2D_flux_basis_r=50", U);

  MatrixDense phi_0(2*mesh->number_cells(), 100);
  Vector k(100);

  // ROM
  double** data = get_samples(100, 13);
  for (int i=0; i<100; i++)
  {
    Material::SP_material mat = get_mat(data[i]);

    ROMSolver<_2D> ROM(inp, mesh, mat);
    SP_vector  ROM_flux;
    ROM_flux = new callow::Vector(2*n, 0.0);
    ROM.Solve(U, ROM_flux);
    double keff_rom = ROM.keff();

    double* flux = &(*ROM_flux)[0];
    phi_0.insert_col(i, flux);
    k[i] = keff_rom;

    std::cout << "ROM k_eff= " << keff_rom << "\n";
   }

  phi_0.print_matlab("IAEA_ROM_flux.txt");
  k.print_matlab("IAEA_ROM_eigenvalue.txt");

 return 0;
}


int test_DEIM_snapshots(int argc, char *argv[])

{
  typedef DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
  typedef detran::DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
  typedef callow::MatrixDense::SP_matrix            SP_matrix;

  InputDB::SP_input inp = get_input();

  inp->get<std::string>("equation") = "diffusion";
  Material::SP_material mat = get_mat(false);
  Mesh2D::SP_mesh mesh = get_mesh();

  SP_lossoperator A (new DiffusionLossOperator(inp, mat, mesh, false, 0.0, false, 1.0, false));
  SP_gainoperator B (new DiffusionGainOperator(inp, mat, mesh, false));

  callow::MatrixDense::SP_matrix  L_snapshots;
  L_snapshots =  new callow::MatrixDense(A->number_nonzeros(), 100);

  callow::MatrixDense::SP_matrix  F_snapshots;
  F_snapshots =  new callow::MatrixDense(B->number_nonzeros(), 100);

  double** data = get_samples(100, 13);
  for (int i=0; i< 100; i++)
  {
    Material::SP_material mat = get_mat(data[i]);

    SP_lossoperator A (new DiffusionLossOperator(inp, mat, mesh, false, 0.0, false, 1.0, false));
    SP_gainoperator B (new DiffusionGainOperator(inp, mat, mesh, false));

    L_snapshots->insert_col(i, A->values());

    F_snapshots->insert_col(i, B->values());
 }

  L_snapshots->print_matlab("L_snapshots.txt");
  F_snapshots->print_matlab("F_snapshots.txt");

return 0;
}

int test_IAEA_DEIM(int argc, char *argv[])

{
  typedef DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
  typedef detran::DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
  typedef callow::MatrixDense::SP_matrix            SP_matrix;
  typedef std::vector<SP_matrix>                    vec_matrix;
  typedef LinearSolverCreator::SP_solver            SP_solver;
  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef callow::EigenSolver::SP_solver            SP_eigensolver;
  typedef callow::EigenSolverCreator                Creator_T;

  double xs[] = {0.030120, 0.080032,  0.130032, 0.02, 0.04, 1.5,
		         0.4, 2, 0.3,  0.135, 0.040160, 0.010024, 0.085032};


  InputDB::SP_input inp = get_input();
  inp->get<std::string>("equation") = "diffusion";
  Material::SP_material mat = get_mat(xs);
  Mesh2D::SP_mesh mesh = get_mesh();

 // 1- offline
  int n = mesh->number_cells();
  SP_lossoperator A (new DiffusionLossOperator(inp, mat, mesh, false, 0.0, false, 1.0, false));
  SP_gainoperator B (new DiffusionGainOperator(inp, mat, mesh, false));

  int rL = 15;
  int rF = 5;
  int r = 50;

  SP_matrix UL;
  UL = new callow::MatrixDense(A->number_nonzeros(), rL);
  ROMBasis::GetBasis("../../../source/src/solvers/test/rom_basis/IAEA2D_L_basis", UL);

  SP_matrix UF;
  UF = new callow::MatrixDense(520000, rF);
  ROMBasis::GetBasis("../../../source/src/solvers/test/rom_basis/IAEA2D_F_basis", UF);

  SP_matrix U;
  U = new callow::MatrixDense(2*n, 2*r);
  ROMBasis::GetBasis("../../../source/src/solvers/test/rom_basis/IAEA2D_flux_basis_r=50", U);

  offline_stage O(A, U, UL, rL);
  vec_matrix M_L;
  M_L = O.Decompositon();
  SP_matrix Udeim_L;
  Udeim_L = O.Ur_deim;
  int* l_L = O.Interpolation_indices();

  offline_stage O_F(B, U, UF, rF);
  vec_matrix M_F;
  M_F = O_F.Decompositon();
  SP_matrix Udeim_F;
  Udeim_F = O_F.Ur_deim;
  int* l_F = O_F.Interpolation_indices();

  SP_solver solver_L;
  solver_L = LinearSolverCreator::Create(inp->template get<SP_input>("inner_solver_db"));
  solver_L->set_operators(Udeim_L, inp->template get<SP_input>("inner_solver_db"));

  SP_solver solver_F;
  solver_F = LinearSolverCreator::Create(inp->template get<SP_input>("inner_solver_db"));
  solver_F->set_operators(Udeim_F, inp->template get<SP_input>("inner_solver_db"));

  SP_eigensolver eigensolver;
  eigensolver = Creator_T::Create(inp->template get<SP_input>("eigen_solver_db"));
  SP_vector  ROM_flux;
  ROM_flux = new callow::Vector(2*n, 0.0);

  MatrixDense phi(2*mesh->number_cells(), 100);
  Vector k_eff(100);

  double** data = get_samples(100, 13);

  for (int sample=0; sample < 100; sample++)
  {
	std::cout << "sample =  " << sample << "\n";

    std::cout << "************** DEIM ****************\n";

    Material::SP_material mat = get_mat(data[sample]);
    SP_lossoperator A (new DiffusionLossOperator(inp, mat, mesh, false, 0.0, false, 1.0, false));
    SP_gainoperator B (new DiffusionGainOperator(inp, mat, mesh, false));

    double* vectorized_L = A->values();
    double* vectorized_F = B->values();

    double* vectorized_F_nnz;
    vectorized_F_nnz= new double[UF->number_rows()];

    int h = 0;
    for (int i=0; i< B->number_nonzeros(); i++)
    {
      if (vectorized_F[i] !=0)
      {
        vectorized_F_nnz[h] = vectorized_F[i];
    	h += 1;
      }
    }

    SP_vector b_L;
    b_L = new callow::Vector(rL, 0.0);

    SP_vector x_L;
	x_L = new callow::Vector(rL, 0.0);

    SP_vector b_F;
    b_F = new callow::Vector(rF, 0.0);

    SP_vector x_F;
	x_F = new callow::Vector(rF, 0.0);

    for (int i=0; i<rL; i++)

    {
      int id = l_L[i];
      (*b_L)(i) = vectorized_L[id];
    }

    solver_L->solve(*b_L, *x_L);

    for (int i=0; i<rF; i++)
    {
      int id = l_F[i];
	  (*b_F)(i) = vectorized_F_nnz[id];
    }

    solver_F->solve(*b_F, *x_F);

    SP_matrix Lr;
    Lr = new MatrixDense(2*r, 2*r);

    double v = 0.0;

    for (int k=0; k< rL; k++)
    {
      for (int i=0; i<2*r; i++)
      {
        for (int j=0; j<2*r; j++)
    	{
    	  v = (*x_L)[k]*(*M_L[k])(i, j);
    	  Lr->insert(i, j, v, 1);
      }
     }
    }

    SP_matrix Fr;
    Fr = new MatrixDense(2*r, 2*r);
    // construct the Left operator
    for (int k=0; k< rF; k++)
    {
      for (int i=0; i<2*r; i++)
      {
        for (int j=0; j<2*r; j++)
        {
          v = (*x_F)[k]*(*M_F[k])(i, j);
          Fr->insert(i, j, v, 1);
        }
      }
    }

    // solve the eigenvalue problem

    eigensolver->set_operators(Fr, Lr);

    SP_vector x_rom;
    SP_vector x;
    x_rom = new Vector(2*r, 0.0);
    x = new Vector(2*r, 1.0);
    eigensolver->solve(x_rom, x);

    double keff =  eigensolver->eigenvalue();

    k_eff[sample] = keff;

    // reconstruct
    int d_n = U->number_rows();

    U->multiply(*x_rom, *ROM_flux);

    double* a = &(*ROM_flux)[0];

    phi.insert_col(sample, a);
  }

  phi.print_matlab("IAEA_DEIM_flux.txt");
  k_eff.print_matlab("IAEA_DEIM_eigenvalue.txt");

  return 0;

}




