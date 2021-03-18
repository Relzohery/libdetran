/*
 * test_DEIM.cc
 *
 *  Created on: Jan 1, 2021
 *      Author: rabab
 */




#define TEST_LIST\
        FUNC(test_DEIM)\
        FUNC(test_offline)

#include "TestDriver.hh"
#include "solvers/rom/ROMBasis.hh"
#include "solvers/rom/ROMSolver.hh"
#include "solvers/rom/DEIM.hh"
#include "utilities/MathUtilities.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/rom/ParametricRom.hh"
#include "projection_fixture.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "callow/solver/EigenSolverCreator.hh"


using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;

typedef callow::MatrixDense::SP_matrix        SP_matrix;
typedef callow::Matrix::SP_matrix             SP_matrix_sparse;
typedef detran::DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
typedef detran::DiffusionGainOperator::SP_gainoperator    SP_gainoperator;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}


int test_DEIM(int argc, char *argv[])
{

 SP_matrix U;
 U = new callow::MatrixDense(976, 976);
 const char* fname= "/home/rabab/Desktop/DEIM_basis_R";
 ROMBasis::GetBasis(fname, U);

 DEIM D(U, 20);
 D.Search();
 int* l = D.interpolation_indices();
 SP_matrix Um;
 Um = new callow::MatrixDense(20, 20);
 Um = D.ReducedBasis();

 DEIM D2(U, 20);
 D2.Search();
 //int* l2 = D2.interpolation_indices();
return 0;
}

int test_offline(int argc, char *argv[])
{
  SP_matrix UL;
  UL = new callow::MatrixDense(976, 976);
  const char* fname1= "/home/rabab/Desktop/DEIM_basis_L";
  ROMBasis::GetBasis(fname1, UL);

  SP_matrix UR;
  UR = new callow::MatrixDense(560, 560);
  const char* fname2= "/home/rabab/Desktop/DEIM_basis_R";
  ROMBasis::GetBasis(fname2, UR);

  SP_matrix Uf;
  Uf = new callow::MatrixDense(280, 14);
  const char* fname3= "/home/rabab/Desktop/flux_basis";
  ROMBasis::GetBasis(fname3, Uf);

  Mesh1D::SP_mesh mesh = get_mesh(1, "core");
  Material::SP_material mat = get_mat();
  InputDB::SP_input input = get_input();

  input->put<std::string>("operator", "diffusion");
  input->put<std::string>("equation", "diffusion");


  ParametricRom<_1D> R(Uf, UL,  UR, input,  mesh, mat, 15);
  SP_matrix_sparse d_A;
  SP_lossoperator A (new DiffusionLossOperator(input, mat, mesh, false, 0.0, false, 1.0));
  SP_gainoperator B (new DiffusionGainOperator(input, mat, mesh, false));

  d_A = A;
  offline_stage off(A, Uf, UL, 2);
  R.solve();

  //off.VectorToMatrix();
  EigenvalueManager<_1D> manager(input, mat, mesh);

  manager.solve();
  double keff_fom = manager.state()->eigenvalue();

  std::cout << keff_fom << "\n";


return 0;
}
