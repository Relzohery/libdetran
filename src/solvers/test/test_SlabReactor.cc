#define TEST_LIST               \
	FUNC(test_SlabReactor_fom)\
   	FUNC(test_Diffusion_Projection)\
	FUNC(test_projection)\
	FUNC(test_EnergyDependent_projection)\


#include "Mesh1D.hh"
#include "Material.hh"
#include "TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/solver/EigenSolverCreator.hh"
#include <vector>
#include <iostream>
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "boundary/BoundarySN.hh"
#include "kinetics/LinearExternalSource.hh"
#include "kinetics/LinearMaterial.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "Eigensolver.hh"

#include "EigenArnoldi.hh"
#include "callow/solver/EigenSolverCreator.hh"
#include "utilities/InputDB.hh"
#include <string>
#include "FixedSourceManager.hh"
#include "solvers/eigen/EnergyDependentEigenLHS.hh"
#include <iostream>
#include <fstream>




using namespace detran_test;
using namespace detran;
using namespace detran_material;
using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;



int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

Material::SP_material get_mat()
{
  Material::SP_material mat = Material::Create(4, 2, "slabreactor");
  // Material 0: Water
  // Total
  mat->set_sigma_t(0, 0, 0.1890);
  mat->set_sigma_t(0, 1, 1.4633);

  // Fission
  mat->set_sigma_f(0, 0, 0.0);
  mat->set_sigma_f(0, 1, 0.0);
  mat->set_chi(0, 0, 0.0);
  mat->set_chi(0, 1, 0.0);

  // Scattering
  mat->set_sigma_s(0, 0, 0, 0.1507);
  mat->set_sigma_s(0, 0, 1, 0.0000);
  mat->set_sigma_s(0, 1, 0, 0.0380);
  mat->set_sigma_s(0, 1, 1, 1.4536);
  mat->compute_diff_coef();                   ;

  // Material 1: Fuel I
  // Total
  mat->set_sigma_t(1, 0, 0.2263);
  mat->set_sigma_t(1, 1, 1.0119);

  // Fission
  mat->set_sigma_f(1, 0, 0.0067);
  mat->set_sigma_f(1, 1, 0.1241);
  mat->set_chi(1, 0, 1.0);
  mat->set_chi(1, 1, 0.0);

  // Scattering
  mat->set_sigma_s(1, 0, 0, 0.2006);
  mat->set_sigma_s(1, 0, 1, 0.0000);
  mat->set_sigma_s(1, 1, 0, 0.0161);
  mat->set_sigma_s(1, 1, 1, 0.9355);
  mat->compute_diff_coef();                   ;


  // Material 3: Fuel II

  // Total
  mat->set_sigma_t(2, 0, 0.2252);
  mat->set_sigma_t(2, 1, 0.9915);

  // Fission
  mat->set_sigma_f(2, 0, 0.0078);
  mat->set_sigma_f(2, 1, 0.1542);
  mat->set_chi(2, 0, 1.0);
  mat->set_chi(2, 1, 0.0);

  // Scattering
  mat->set_sigma_s(2, 0, 0, 0.1995);
  mat->set_sigma_s(2, 0, 1, 0.0000);
  mat->set_sigma_s(2, 1, 0, 0.0156);
  mat->set_sigma_s(2, 1, 1, 0.9014);
  mat->compute_diff_coef();                   ;


  // Material 4: Fuel II + Gd

 // Total
  mat->set_sigma_t(3, 0, 0.2173);
  mat->set_sigma_t(3, 1, 1.0606);

  // Fission
  mat->set_sigma_f(3, 0, 0.0056);
  mat->set_sigma_f(3, 1, 0.0187);
  mat->set_chi(3, 0, 1.0);
  mat->set_chi(3, 1, 0.0);

   // Scattering
  mat->set_sigma_s(3, 0, 0, 0.1902);
  mat->set_sigma_s(3, 0, 1, 0.0000);
  mat->set_sigma_s(3, 1, 0, 0.0136);
  mat->set_sigma_s(3, 1, 1, 0.5733);
  mat->compute_diff_coef();                   ;



  mat->finalize();

  return mat;
}

Mesh1D::SP_mesh get_mesh(Mesh1D::size_t fmm = 1)
{
  vec_dbl cm_assembly(7, 0.0);
  cm_assembly[0] = 0.0;
  cm_assembly[1] = 1.1580;
  cm_assembly[2] = 4.4790;
  cm_assembly[3] = 7.8000;
  cm_assembly[4] = 11.1210;
  cm_assembly[5] = 14.4420;
  cm_assembly[6] = 15.6000;

  std::vector<double> cm_core;
  std::vector<double> fm_core;

  std::vector<int> fm_assembly = {2, 4, 4, 4, 4, 2};
  for (int i=0; i < 7;  i++)
  {
    for (int j=0; j <7; j++)
    {
     cm_core.push_back(cm_assembly[j] + 15.6*i);
     fm_core.push_back(fm_assembly[j]);
    }
  }

  cm_core.assign(0.0 ,0.0);


  Mesh1D::vec_int mt_0(6);
  mt_0[0] = 0; mt_0[1] =1;
  mt_0[2] = 2; mt_0[3] = 2;
  mt_0[4] = 1; mt_0[5] = 0;

  Mesh1D::vec_int mt_1(6);
  mt_1[0] = 0; mt_1[1] =1;
  mt_1[2] = 1; mt_1[3] = 1;
  mt_1[4] = 1; mt_1[5] = 0;

  Mesh1D::vec_int mt_2(6);
  mt_2[0] = 0; mt_2[1] =1;
  mt_2[2] = 3; mt_2[3] = 3;
  mt_2[4] = 1; mt_2[5] = 0;

  Mesh1D::vec_int mt_3(6);
  mt_3[0] = 0; mt_3[1] = 3;
  mt_3[2] = 3; mt_3[3] = 3;
  mt_3[4] = 3; mt_3[5] = 0;

  Mesh1D::SP_mesh mesh = Mesh1D::Create(fm_assembly, cm_assembly, mt_0);

  return mesh;

}

InputDB::SP_input get_input()
{
  InputDB::SP_input inp(new InputDB("Slab Reactor"));
  inp->put<std::string>("problem_type", "eigenvalue");
  inp->put<int>("number_groups",                  2);
  inp->put<std::string>("equation",               "dd");
  inp->put<std::string>("inner_solver",           "SI");
  inp->put<int>("inner_max_iters",            1000);
  inp->put<int>("inner_max_iters",            1000);
  inp->put<double>("inner_tolerance",            1e-7);
  inp->put<int>("inner_print_level",          0);
  inp->put<std::string>("outer_solver",              "GMRES");
  inp->put<int>("outer_max_iters",            1000);
  inp->put<double>("outer_tolerance",            1e-7);
  inp->put<int>("outer_print_level",          0);
  inp->put<std::string>("eigen_solver",       "arnoldi");
  inp->put<int>("eigen_max_iters",            200);
  inp->put<double>("eigen_tolerance",            1e-7);
  inp->put<std::string>("bc_west",                    "reflect");
  inp->put<std::string>("bc_east",                    "reflect");
  inp->put<int>("quad_number_polar_octant",   16);

  return inp;
}

int test_SlabReactor_fom(int argc, char *argv[])
{
  InputDB::SP_input inp = get_input();
  Mesh1D::SP_mesh mesh = get_mesh();
  Material::SP_material mat = get_mat();
  EigenvalueManager<_1D> manager(inp, mat, mesh);
  manager.solve();
  State::SP_state ic = manager.state();
  std::cout << "manager.state()->eigenvalue()" << "\n";

  TEST(soft_equiv(1.329576914, manager.state()->eigenvalue(), 1.0e-4));

  return 0;
}

int test_projection(int argc, char *argv[])
{
  typedef typename Eigensolver<_1D>::Fixed_T              Fixed_T;
  typedef typename Eigensolver<_1D>::SP_solver            SP_solver;
  typedef typename Fixed_T::SP_manager                    SP_mg_solver;
  typedef EnergyIndependentEigenOperator<_1D>             Operator_T;
  typedef typename Operator_T::SP_operator                SP_operator;
  typedef EnergyDependentEigenLHS<_1D>                    LHS_Operator_T;
  typename LHS_Operator_T::SP_operator B;


  InputDB::SP_input input = get_input();
  input->put<std::string>("outer_solver", "GMRES");
  input->put<std::string>("eigen_solver", "GD");

  Mesh1D::SP_mesh mesh = get_mesh();
  Material::SP_material mat = get_mat();

  SP_mg_solver mg_solver;
  mg_solver = new FixedSourceManager<_1D>(input, mat, mesh, false, true);
  mg_solver->setup();
  mg_solver->set_solver();

  SP_operator A;
  A = new Operator_T(mg_solver);
  A->compute_explicit("EnergyIndependentEigenOperator");

  int m = A->number_columns();
  int n = A->number_rows();

  // get the basis
  ifstream infile;
  infile.open("/home/rabab/Research/detran_demo/steady_state/"
              "slab_reactor/fission_density_basis_assem0", ios::binary | ios::in);
  int r = 5;

  double U[20][r];

  infile.seekg(0);
  infile.read((char *) &U, sizeof(U)); // read the number of element

  callow::MatrixDense A_r1(n, 5);

  vector<vector<double>> basis_vecs;

  for (int i=0; i<r; i++)
  {
   vector<double> v;
   std::cout << i << "\n";
   for (int j=0; j< n; j++)
   {
     std::cout << U[j][i] << "\n";
     v.push_back(U[j][i]);
   }
    basis_vecs.push_back(v);
  }

  callow::Vector y(n, 0.0);
  for (int i=0; i<r; i++)
  {
    std::cout << i << "\n";
    A->multiply(basis_vecs[i], y);
    A_r1.insert_col(i, y, 0);
  }

  A_r1.assemble();
  A_r1.print_matlab("reduced_Ar11");

 // now do the multiplication by UT
  callow::MatrixDense::SP_matrix A_r;
  A_r = new callow::MatrixDense(5, 5);
  for (int i=0; i<5; i++)
  {
    for (int k=0; k<5 ; k++)
    {
     for (int j=0; j<20; j++)
     {
       double v = U[j][i]*A_r1(j, k);
       //std::cout << U[j][i] << " & ";
       A_r->insert(i, k, v, 1);
     }
    }
  }
 A_r->print_matlab("reduced_A");

return 0;

}

int test_EnergyDependent_projection(int argc, char *argv[])
{
  typedef typename Eigensolver<_1D>::Fixed_T              Fixed_T;
  typedef typename Eigensolver<_1D>::SP_solver            SP_solver;
  typedef typename Fixed_T::SP_manager                  SP_mg_solver;
  typedef EnergyDependentEigenLHS<_1D>                LHS_Operator_T;
  typename LHS_Operator_T::SP_operator B;

  InputDB::SP_input input = get_input();
  input->put<std::string>("outer_solver", "GMRES");
  input->put<std::string>("eigen_solver", "GD");

  Mesh1D::SP_mesh mesh = get_mesh();
  Material::SP_material mat = get_mat();

  SP_mg_solver mg_solver_ED;
  mg_solver_ED = new FixedSourceManager<_1D>(input, mat, mesh, false, true);
  mg_solver_ED->setup();
  mg_solver_ED->set_solver();

  B = new LHS_Operator_T(mg_solver_ED);
  int m = B->number_columns();
  int n = B->number_rows();

  std::cout << m << " value of m" << "\n";
  std::cout << n << " value of n" << "\n";

  // get the basis
  ifstream infile;
  infile.open("/home/rabab/Research/detran_demo/steady_state/"
              "slab_reactor/flux_basis_assem0_diff", ios::binary | ios::in);
  int r = 10;

  double U[n][r];

  infile.seekg(0);
  infile.read((char *) &U, sizeof(U)); // read the number of element

  callow::MatrixDense B_r1(n, r);

  vector<vector<double>> basis_vecs;

  for (int i=0; i<r; i++)
  {
   vector<double> v;
   for (int j=0; j< n; j++)
   {
     std::cout << U[j][i] << "\n";
     v.push_back(U[j][i]);
   }
    basis_vecs.push_back(v);
  }

  callow::Vector y(n, 0.0);
  for (int i=0; i<r; i++)
  {
    B->multiply(basis_vecs[i], y);
    B_r1.insert_col(i, y, 0);
  }

 // now do the multiplication by UT
  callow::MatrixDense::SP_matrix B_r;
  B_r = new callow::MatrixDense(r, r);
  for (int i=0; i<r; i++)
  {
    for (int k=0; k<r ; k++)
    {
     for (int j=0; j<m; j++)
     {
       double v = U[j][i]*B_r1(j, k);
       B_r->insert(i, k, v, 1);
     }
    }
  }
  B_r->print_matlab("reduced_B");

return 0;
}

int test_Diffusion_Projection(int argc, char *argv[])
{
 typedef DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
 typedef DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
 typedef callow::EigenSolver::SP_solver            SP_eigensolver;
 typedef callow::EigenSolverCreator                Creator_T;
 typedef callow::Vector::SP_vector                 SP_vector;

 InputDB::SP_input inp = get_input();
 Mesh1D::SP_mesh mesh = get_mesh();
 Material::SP_material mat = get_mat();
 inp->put<std::string>("equation",  "diffusion");

 EigenvalueManager<_1D> manager(inp, mat, mesh);
 manager.solve();
 State::SP_state ic = manager.state();

 std::cout << "manager.state()->eigenvalue()" << "\n";

 SP_lossoperator A (new DiffusionLossOperator(inp, mat, mesh, false, 0.0, false, 1.0));
 SP_gainoperator B (new DiffusionGainOperator(inp, mat, mesh, false));

 A->print_matlab("loss_operator");
 B->print_matlab("gain_operator");

 // get the basis
 ifstream infile;

 infile.open("/home/rabab/Research/detran_demo/steady_state/"
               "slab_reactor/flux_basis_assem0_diff", ios::binary | ios::in);

 int m = A->number_columns();
 int n = A->number_rows();

 int r = 10;
 double U[n][r];

 infile.seekg(0);
 infile.read((char *) &U, sizeof(U)); // read the number of element

 callow::MatrixDense U_mat(n, r);

 callow::MatrixDense A_r1(n, 10);
 callow::MatrixDense B_r1(n, 10);

 vector<vector<double>> basis_vecs;

 for (int i=0; i<r; i++)
 {
   vector<double> v;
   for (int j=0; j< n; j++)
   {
     v.push_back(U[j][i]);
   }
     basis_vecs.push_back(v);
 }

 callow::Vector y(n, 0.0);
 for (int i=0; i<r; i++)
 {
   A->multiply(basis_vecs[i], y);
   A_r1.insert_col(i, y, 0);
   U_mat.insert_col(i, basis_vecs[i], 0);

   B->multiply(basis_vecs[i], y);
   B_r1.insert_col(i, y, 0);
 }

 // now do the multiplication by UT
 callow::MatrixDense::SP_matrix A_r;
 callow::MatrixDense::SP_matrix B_r;

 A_r = new callow::MatrixDense(10, 10);
 B_r = new callow::MatrixDense(10, 10);

 for (int i=0; i<10; i++)
 {
   for (int k=0; k<10 ; k++)
   {
     for (int j=0; j<40; j++)
     {
       double va = U[j][i]*A_r1(j, k);
       double vb = U[j][i]*B_r1(j, k);
       A_r->insert(i, k, va, 1);
       B_r->insert(i, k, vb, 1);
     }
   }
 }

  A_r->print_matlab("reduced_diff_A");
  B_r->print_matlab("reduced_diff_B");

  SP_eigensolver eigensolver;
  eigensolver = Creator_T::Create(inp);
  eigensolver->set_operators(B_r, A_r, inp);

  SP_vector x;
  SP_vector x1;
  x  = new callow::Vector(10, 0.0);
  x1 = new callow::Vector(10, 1.0);
  eigensolver->solve(x, x1);
  std::cout << eigensolver->eigenvalue() << "\n";

  // reconstruct flux
  callow::Vector phi0(20, 0.0);
  callow::Vector phi1(20, 0.0);
  callow::Vector phi(m, 0.0);

  std::cout <<"the norm is  " << x->norm() << "\n";
  double norm = x->norm();

  U_mat.multiply(*x, phi);
  for (int i=0; i<r ; i++)
  {
    std::cout << (*x)[i]/norm << "\n";
  }

  for (int i=0; i<20; i++)
  {
   phi0[i] = phi[i];
   phi1[i] = phi[i+10];
  }

  return 0;
}
