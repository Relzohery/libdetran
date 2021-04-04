/*
 * test_IAEA2D.cc
 *
 *  Created on: Mar 22, 2021
 *      Author: rabab
 */


#define TEST_LIST \
        FUNC(test_IAEA2D)\
		FUNC(test_IAEA2D_ROM)\
		FUNC(test_DEIM_snapshots)


#include "TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "utilities/MathUtilities.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "solvers/EigenvalueManager.hh"
#include "callow/vector/Vector.hh"
#include "solvers/rom/ROMSolver.hh"
#include "Mesh2D.hh"
#include "DEIM.hh"
#include "offline_stage.hh"
#include <random>

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


//typedef callow::MatrixDense::SP_matrix SP_matrix;
//typedef callow::Vector::SP_vector      SP_vector;



int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}


InputDB::SP_input get_input()
 {
   InputDB::SP_input inp(new InputDB("IAEA-2D"));
   inp->put<int>("number_groups",                      2);
   inp->put<int>("dimension",                          2);
   inp->put<std::string>("equation",                           "diffusion");
   inp->put<std::string>("bc_west",                            "reflect");
   inp->put<std::string>("bc_east",                            "vacuum");
   inp->put<std::string>("bc_south",                           "reflect");
   inp->put<std::string>("bc_north",                           "vacuum");
   inp->put<int>("eigen_max_iters",                    1000);
   inp->put<std::string>("eigen_solver",                       "arnoldi");
   InputDB::SP_input db(new InputDB("callow dp"));
   db->put<double>("linear_solver_atol",              1e-9);
   db->put<double>("linear_solver_rtol",                  1e-8);
   db->put<std::string>("linear_solver_type",                  "petsc");
   db->put<int>("linear_solver_maxit",                 5000);
   db->put<int>("linear_solver_gmres_restart",         30);
   db->put<std::string>("eigen_solver_type",                   "slepc");
   db->put<int>("linear_solver_monitor_level",         0);
   db->put<InputDB::SP_input>("outer_solver_db", db);
   db->put<InputDB::SP_input>("eigen_solver_db", db);

   return inp;

 }

Material::SP_material get_mat(bool perturb = false)
{


 std::random_device rd;  //Will be used to obtain a seed for the random number engine
 std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
 std::uniform_real_distribution<double> distribution(0.0,1.0);

 double number = distribution(generator)/100;

  double xs[] ={0.030120, 0.080032,  0.130032, 0.02, 0.04,  1.5, 0.4, 2, 0.3,  0.135, 0.040160, 0.010024, 0.085032};

  if (perturb)
  {
	for (int i=0; i<13; i++)
	{
	  number = distribution(generator)/100;
	  xs[i] *= number;
	}
  }


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
 mat->set_sigma_t(1, 0, xs[2]);
 mat->set_sigma_t(1, 1, xs[12]);
 mat->set_sigma_s(1, 1, 0, xs[3]);
 mat->set_sigma_f(1, 1, xs[9]);
 mat->set_diff_coef(1, 0, xs[5]);
 mat->set_diff_coef(1, 1, xs[6]);
 mat->set_chi(1, 0, 1.0);
 mat->compute_sigma_a();

 // Material 2
 mat->set_sigma_t(2, 0, xs[0]);
 mat->set_sigma_t(2, 1, xs[1]);
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
 mat->display();

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

int test_IAEA2D(int argc, char *argv[])
{
  time_t begin, end;
  time(&begin);

  InputDB::SP_input inp = get_input();

  //-------------------------------------------------------------------------//
  // MATERIAL
  //-------------------------------------------------------------------------//
  inp->get<std::string>("equation") = "diffusion";

  Material::SP_material mat = get_mat(false);

  Mesh2D::SP_mesh mesh = get_mesh();

  EigenvalueManager<_2D> manager(inp, mat, mesh);
  manager.solve();
  State::SP_state ic = manager.state();
}


int test_IAEA2D_ROM(int argc, char *argv[])
{
  time_t begin, end;
  time(&begin);
  std::cout << " input **********\n";
  InputDB::SP_input inp = get_input();

  //-------------------------------------------------------------------------//
  // MATERIAL
  //-------------------------------------------------------------------------//
  inp->get<std::string>("equation") = "diffusion";

  Material::SP_material mat = get_mat();

  Mesh2D::SP_mesh mesh = get_mesh();

  EigenvalueManager<_2D> manager(inp, mat, mesh);
  manager.solve();
  State::SP_state ic = manager.state();
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

  double* nnzeros = A->values();
  callow::MatrixDense::SP_matrix  L_snapshots;
  L_snapshots =  new callow::MatrixDense(A->number_nonzeros(), 100);

  callow::MatrixDense::SP_matrix  F_snapshots;
  F_snapshots =  new callow::MatrixDense(B->number_nonzeros(), 100);

  for (int i=0; i< 100; i++)
  {
    Material::SP_material mat = get_mat(true);

    SP_lossoperator A (new DiffusionLossOperator(inp, mat, mesh, false, 0.0, false, 1.0, false));
    SP_gainoperator B (new DiffusionGainOperator(inp, mat, mesh, false));

    L_snapshots->insert_col(i, A->values());

    F_snapshots->insert_col(i, B->values());

	//std::cout << d_M->number_nonzeros() << "\n";
	//LossMatrix_snaps = new callow::MatrixDense(d_M->number_nonzeros(), d_number_steps+1);
 }

  L_snapshots->print_matlab("L_snapshots");
  F_snapshots->print_matlab("F_snapshots");

return 0;
}




