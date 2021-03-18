//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_LRA.cc
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  Implements the LRA kinetics benchmark
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST               \
        FUNC(test_LRA)        \
		FUNC(test_LRA_ROM)

#include "TestDriver.hh"
#include "TimeStepper.hh"
#include "Mesh2D.hh"
#include "external_source/ConstantSource.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "boundary/BoundaryTraits.hh"
#include "boundary/BoundarySN.hh"
#include "kinetics/LinearExternalSource.hh"
#include "kinetics/LinearMaterial.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/time/LRA.hh"
#include "TransientSolver.hh"
#include "ROMBasis.hh"

//
//#include "angle/test/quadrature_fixture.hh"
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
//#include "external_source/test/external_source_fixture.hh"

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

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void test_monitor(void* data, TimeStepper<_2D>* ts, int step, double t,
                  double dt, int it, bool conv)
{
  static double maxp = 0;
  double F = 0;
  TimeStepper<_2D>::vec_int matmap = ts->mesh()->mesh_map("MATERIAL");
  TimeStepper<_2D>::vec_int cmm = ts->mesh()->mesh_map("COARSEMESH");

  vec_dbl &T = ts->multiphysics()->variable(0);
  //TimeStepper<_2D>::size_t N = 0;
  double T_avg = 0;
  double T_max = 0;
  double V = 0;
  for (int i = 0; i < ts->mesh()->number_cells(); ++i)
  {
    //printf(" %5i %18.12e \n", i, T[i]);
    int m = matmap[i];
    F += ts->mesh()->volume(0) *
         (ts->state()->phi(0)[i] * ts->material()->sigma_f(m, 0) +
          ts->state()->phi(1)[i] * ts->material()->sigma_f(m, 1) );
    if (cmm[i] != 4)
    {
      V += ts->mesh()->volume(0);
      T_avg += ts->mesh()->volume(0) * T[i];
      if (T[i] > T_max) T_max = T[i];
    }
  }
  F *= detran_user::KAPPA / V;
  if (F > maxp && conv) maxp = F;
  T_avg /= 17550.0;
//  std::cout << " phi0=" << ts->state()->phi(0)[0] << " "
//            << " phi1=" << ts->state()->phi(1)[0] << std::endl;
  printf("** %5i  %16.13f  %18.12e  %18.12e  %18.12e  %5i \n", step, t, F, T_avg, T_max, it);
}

//---------------------------------------------------------------------------//
Mesh2D::SP_mesh get_mesh(Mesh2D::size_t fmm = 1)
{
  Mesh2D::vec_dbl cm(12, 0.0);
  for (int i = 1; i < 12; ++i) cm[i] = cm[i-1] + 15.0;
  Mesh2D::vec_int fm(11, fmm);

  int icmm[]={1 , 0 , 0 , 0 , 0 , 1 , 1 , 2 , 2 , 4 , 4,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 4 , 4,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 4 , 4,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 4 , 4,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 4 , 4,
              1 , 0 , 0 , 0 , 0 , 1 , 1 , 5 , 5 , 4 , 4,
              1 , 0 , 0 , 0 , 0 , 1 , 1 , 5 , 5 , 4 , 4,
              2 , 2 , 2 , 2 , 2 , 2 , 2 , 3 , 4 , 4 , 4,
              2 , 2 , 2 , 2 , 2 , 2 , 2 , 4 , 4 , 4 , 4,
              4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4,
              4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4};
  Mesh2D::vec_int cmm(121);
  for (int i = 0; i < 121; ++i) cmm[i] = icmm[i];

  int ibmm[]={0,  1,  2,  3,  4,  5,  6,  7,  8, 78, 78,
              9, 10, 11, 12, 13, 14, 15, 16, 17, 78, 78,
             18, 19, 20, 21, 22, 23, 24, 25, 26, 78, 78,
             27, 28, 29, 30, 31, 32, 33, 34, 35, 78, 78,
             36, 37, 38, 39, 40, 41, 42, 43, 44, 78, 78,
             45, 46, 47, 48, 49, 50, 51, 52, 53, 78, 78,
             54, 55, 56, 57, 58, 59, 60, 61, 62, 78, 78,
             63, 64, 65, 66, 67, 68, 69, 70, 78, 78, 78,
             71, 72, 73, 74, 75, 76, 77, 78, 78, 78, 78,
             78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
             78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78};
  Mesh2D::vec_int bmm(121);
  for (int i = 0; i < 121; ++i) bmm[i] = ibmm[i];

  Mesh2D::SP_mesh mesh = Mesh2D::Create(fm, fm, cm, cm, cmm);
  mesh->add_coarse_mesh_map("ASSEMBLY",    bmm);
  mesh->add_coarse_mesh_map("SUBASSEMBLY", bmm);
  mesh->add_coarse_mesh_map("COARSEMESH",  cmm);

  // fine mesh material map
  Mesh2D::vec_int  mm(mesh->number_cells(), 0);
  for (int i = 0; i < mm.size(); ++i)
    mm[i]  = i;
  mesh->add_mesh_map("MATERIAL", mm);
  return mesh;
}


//---------------------------------------------------------------------------//
InputDB::SP_input get_input()

{     typedef TimeStepper<_2D> TS_2D;
      typedef std::string str;
	  //-------------------------------------------------------------------------//
	  // INPUT
	  //-------------------------------------------------------------------------//


InputDB::SP_input inp(new InputDB("LRA benchmark"));
  inp->put<int>("number_groups",                  2);
  inp->put<int>("dimension",                      2);
  inp->put<str>("equation",                       "diffusion");
  inp->put<str>("bc_west",                        "reflect");
  inp->put<str>("bc_east",                        "vacuum");
  inp->put<str>("bc_south",                       "reflect");
  inp->put<str>("bc_north",                       "vacuum");
  inp->put<int>("bc_zero_flux",                   0);
  inp->put<int>("quad_number_polar_octant",       3);
  inp->put<int>("quad_number_azimuth_octant",     3);
  inp->put<str>("eigen_solver",                   "arnoldi");
  inp->put<double>("eigen_tolerance",                1e-12);
  inp->put<int>("eigen_max_iters",                1000);
  inp->put<str>("outer_solver",                   "GMRES");
  inp->put<double>("outer_tolerance",                1e-12);
  inp->put<int>("outer_max_iters",                1000);
  inp->put<int>("outer_print_level",              0);
  inp->put<int>("outer_krylov_group_cutoff",      0);
  inp->put<str>("outer_pc_type",                  "mgdsa");
  inp->put<str>("inner_solver",                   "GMRES");
  inp->put<double>("inner_tolerance",                1e-12);
  inp->put<int>("inner_max_iters",                1000);
  inp->put<int>("inner_print_level",              0);
  inp->put<str>("inner_pc_type",                  "DSA");
  // gmres parameters
  InputDB::SP_input db(new InputDB("callow_db"));
   db->put<double>("linear_solver_atol",              0.0);
  db->put<double>("linear_solver_rtol",              1e-16);
   db->put<str>("linear_solver_type",              "gmres");
  db->put<string>("petsc_pc_type",                    "lu");
  db->put<int>("linear_solver_maxit",             1000);
  db->put<int>("linear_solver_gmres_restart",     30);
  db->put<int>("linear_solver_monitor_level",     0);
  db->put<str>("pc_type",                         "petsc_pc");
  //db->put<str>("petsc_pc_type",                   "lu");
  db->put<int>("petsc_pc_factor_levels",          3);
  db->put<str>("eigen_solver_type",               "power");
  db->put<int>("eigen_solver_maxit",              1000);
  db->put<int>("eigen_solver_monitor_level",      1);
  db->put<double>("eigen_solver_tol",                1.0e-12);
  inp->put<InputDB::SP_input>("inner_solver_db", db);
  inp->put<InputDB::SP_input>("inner_pc_db", db);
  inp->put<InputDB::SP_input>("outer_solver_db", db);
  inp->put<InputDB::SP_input>("eigen_solver_db", db);
  inp->put<int>("ts_max_steps",                   10000);
  //inp->put<int>("ts_scheme",                      TS_2D::BDF2);
  inp->put<int>("ts_output",                      0);
  inp->put<double>("ts_step_size",                   0.001);
  inp->put<double>("ts_final_time",                  3.0);
  //inp->put<int>("ts_no_extrapolation",            1);
  inp->put<int>("ts_max_iters",                   10);
  inp->put<double>("ts_tolerance",                   1.0e-8);
  #
  InputDB::SP_input preconditioner_db(new InputDB("preconditioner_db"));
  preconditioner_db->put<double>("linear_solver_atol",              0.0);
  preconditioner_db->put<double>("linear_solver_rtol",              1e-12);
  preconditioner_db->put<str>("linear_solver_type",              "gmres");
  preconditioner_db->put<int>("linear_solver_maxit",             5000);
  preconditioner_db->put<int>("linear_solver_gmres_restart",     30);
  preconditioner_db->put<int>("linear_solver_monitor_level",     0);
  preconditioner_db->put<str>("pc_type",                         "petsc_pc");
  preconditioner_db->put<str>("petsc_pc_type",                   "ilu");
  preconditioner_db->put<int>("petsc_pc_factor_levels",          2);


  InputDB::SP_input rom_db(new InputDB("inner_solver_db"));
  rom_db->put<string>("linear_solver_type",                 "petsc");
  rom_db->put<string>("pc_type",                            "petsc_pc");
  rom_db->put<double>("linear_solver_rtol",              1e-16);
  rom_db->put<string>("petsc_pc_type",                      "lu");
  rom_db->put<int>("linear_solver_maxit",                   1000);
  rom_db->put<int>("linear_solver_gmres_restart",           30);
  rom_db->put<int>("linear_solver_monitor_level",           0);
  rom_db->put<string>("eigen_solver_type",                  "slepc");
  rom_db->put<int>("petsc_pc_factor_levels",          3);
  db->put<double>("eigen_solver_tol",                   1e-15);
  inp->put<InputDB::SP_input>("inner_solver_db",        db);
  inp->put<InputDB::SP_input>("outer_solver_db",        db);
  inp->put<InputDB::SP_input>("eigen_solver_db",        db);
  inp->put<InputDB::SP_input>("rom_solver_db",          db);

  return inp;

}

//---------------------------------------------------------------------------//
int test_LRA(int argc, char *argv[])
{

  typedef TimeStepper<_2D> TS_2D;
  typedef callow::MatrixDense::SP_matrix            SP_matrix;
  typedef callow::Vector::SP_vector                 SP_vector;
   InputDB::SP_input inp = get_input();
  //-------------------------------------------------------------------------//
  // MESH
  //-------------------------------------------------------------------------//

  TS_2D::SP_mesh mesh = get_mesh(2);

  //-------------------------------------------------------------------------//
  // MATERIAL
  //-------------------------------------------------------------------------//

  bool transport = false;
  if (inp->get<std::string>("equation") != "diffusion") transport = true;
  TS_2D::SP_material mat(new detran_user::LRA(mesh, transport, false));

  //-------------------------------------------------------------------------//
  // STEADY STATE
  //-------------------------------------------------------------------------//

  EigenvalueManager<_2D> manager(inp, mat, mesh);
  manager.solve();
  State::SP_state ic = manager.state();
  mat->set_eigenvalue(ic->eigenvalue());
  mat->update(0, 0, 1, false);

  // Normalize state.
  double F = 0;
  TS_2D::vec_int matmap = mesh->mesh_map("MATERIAL");
  for (int i = 0; i < mesh->number_cells(); ++i)
  {
    int m = matmap[i];
    F += mesh->volume(0) *
         (ic->phi(0)[i] * mat->sigma_f(m, 0) +
          ic->phi(1)[i] * mat->sigma_f(m, 1));
  }
  F *= detran_user::KAPPA / 17550.0;
  ic->scale(1.0e-6/F);



  //-------------------------------------------------------------------------//
  // TIME STEPPER
  //-------------------------------------------------------------------------//

  TS_2D stepper(inp, mat, mesh, true);
  stepper.set_monitor(test_monitor);
  detran_utilities::SP<detran_user::LRA> mat_lra;
  mat_lra = mat;

  stepper.set_multiphysics(mat_lra->physics(),
                           detran_user::update_T_rhs<_2D>,
                           (void *) mat_lra.bp());

  stepper.solve(ic);



  printf(" %20.16f %20.16f ", ic->phi(0)[0], ic->phi(0)[1]);
  std::cout << std::endl;

  State::SP_state final = stepper.state();
  SP_matrix phi_mat;
  SP_matrix precursors_mat;
  SP_matrix power_mat;
  SP_vector power;

  phi_mat = stepper.flux_mat;
  precursors_mat = stepper.precursors_mat;
  power_mat = stepper.power_mat;
  power = stepper.power;

  phi_mat->print_matlab("lra_flux_fine.txt");
  precursors_mat->print_matlab("lra_precursors_fine.txt");
  power_mat->print_matlab("lra_spatail_power_fine.txt");
  power->print_matlab("lra_power_fine.txt");

  printf(" %20.16f %20.16f ", final->phi(0)[0], final->phi(0)[1]);
  std::cout << std::endl;

  return 0;

}

//---------------------------------------------------------------------------//
int test_LRA_ROM(int argc, char *argv[])
{

  typedef TimeStepper<_2D> TS_2D;
  typedef callow::MatrixDense::SP_matrix            SP_matrix;


   InputDB::SP_input inp = get_input();
  //-------------------------------------------------------------------------//
  // MESH
  //-------------------------------------------------------------------------//

  TS_2D::SP_mesh mesh = get_mesh(2);

  //-------------------------------------------------------------------------//
  // MATERIAL
  //-------------------------------------------------------------------------//

  bool transport = false;
  if (inp->get<std::string>("equation") != "diffusion") transport = true;
  TS_2D::SP_material mat(new detran_user::LRA(mesh, transport, false));

  //-------------------------------------------------------------------------//
  // STEADY STATE
  //-------------------------------------------------------------------------//

  EigenvalueManager<_2D> manager(inp, mat, mesh);
  manager.solve();
  State::SP_state ic = manager.state();
  mat->set_eigenvalue(ic->eigenvalue());
  mat->update(0, 0, 1, false);

  // Normalize state.
  double F = 0;
  TS_2D::vec_int matmap = mesh->mesh_map("MATERIAL");
  for (int i = 0; i < mesh->number_cells(); ++i)
  {
    int m = matmap[i];
    F += mesh->volume(0) *
         (ic->phi(0)[i] * mat->sigma_f(m, 0) +
          ic->phi(1)[i] * mat->sigma_f(m, 1));
  }
  F *= detran_user::KAPPA / 17550.0;
  ic->scale(1.0e-6/F);

  //-------------------------------------------------------------------------//
  // TIME STEPPER
  //-------------------------------------------------------------------------//

  const char* flux_basis = "./../../../source/src/solvers/test/rom_basis/lra_flux_basis";
  const char* precursors_basis = "./../../../source/src/solvers/test/rom_basis/lra_precursors_basis";

  int r = 30;
  int n = 484;

  SP_matrix basis_f;
  basis_f = new callow::MatrixDense(n*2, 2*r);
  ROMBasis::GetBasis(flux_basis, basis_f);

  SP_matrix basis_p;
  basis_p = new callow::MatrixDense(n*2, r);
  ROMBasis::GetBasis(precursors_basis, basis_p);

  TransientSolver R(inp, mesh, mat, basis_f, basis_p, basis_p, false);

  detran_utilities::SP<detran_user::LRA> mat_lra;

  mat_lra = mat;

  R.set_multiphysics(mat_lra->physics(),
                     detran_user::update_T_rhs<_2D>,
                     (void *) mat_lra.bp());


  R.Solve(ic);
  //time(&end);
  //time_t elapsed = end - begin;
  //printf("time elapsed %1.6f\n", elapsed);

  return 0;
}
//---------------------------------------------------------------------------//
//              end of test_LRA.cc
//---------------------------------------------------------------------------//
