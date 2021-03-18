//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TransientSolver.cc
 *  @brief TransientSolver class definition.
 *  @note  Copyright(C) 2020 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#include "TransientSolver.hh"


TransientSolver::TransientSolver(SP_input inp, SP_mesh mesh, SP_material material,
		                 SP_matrix flux_basis, SP_matrix precursors_basis, SP_matrix Udeim, bool deim)
:d_mesh(mesh),
 d_material(material),
 d_inp(inp),
 d_precursors_basis(precursors_basis),
 d_flux_basis(flux_basis),
 d_deim_basis(Udeim),
 deim_flag (deim)
{
  d_num_cells = d_mesh->number_cells();
  d_number_groups = d_material->number_groups();
  d_precursors_group = d_material->number_precursor_groups();

  d_rf = d_flux_basis->number_columns();
  d_rc = d_precursors_basis->number_columns();

  if (d_inp->check("ts_final_time"))
    d_final_time = d_inp->template get<double>("ts_final_time");
  Assert(d_final_time > 0.0);

  if (d_inp->check("ts_step_size"))
    d_dt = d_inp->template get<double>("ts_step_size");
  Assert(d_dt > 0.0);

  // Compute the number of steps.  May result in longer time than requested!
  d_number_steps = std::ceil(d_final_time / d_dt);

  d_P0 = new callow::Vector (d_num_cells* d_precursors_group, 0.0);
  d_phi0 = new callow::Vector(d_num_cells* d_number_groups, 0.0);

  d_phi0_r = new callow::Vector(d_rf, 0.0);
  d_P0_r = new callow::Vector (d_rc, 0.0);

  // matrix of precursors concentraion at all time steps
  d_precursors = new callow::MatrixDense(d_num_cells*d_precursors_group, d_number_steps+1);
  // matrix of the flux at all time steps
  d_flux = new callow::MatrixDense(d_number_groups*d_num_cells, d_number_steps+1);

  d_A = new callow::MatrixDense(d_rf + d_rc, d_rf+d_rc);
  d_A_ = new callow::MatrixDense(d_rf + d_rc, d_rf+d_rc);


  // long vector of flux and precursors
  d_sols = new callow::MatrixDense(d_rf + d_rc, d_number_steps+1);

  // long vector of reduced flux and precursors
  d_sol_r = new callow::Vector(d_rf+d_rc, 0.0);

  d_sol0_r = new callow::Vector(d_rf + d_rc, 0.0);

  if (d_inp->check("rom_solver_db"))
  {
    db = d_inp->template get<SP_input>("rom_solver_db");

  }
  else db = inp;

  d_solver = LinearSolverCreator::Create(db);
}

//------------------------------------------------------------------------------------//
void TransientSolver::initialize_precursors()
{
  d_fissionsource = new FissionSource(d_state, d_mesh, d_material);
  d_fissionsource->update();

  const State::moments_type &fd = d_fissionsource->density();

  const vec_int &mt = d_mesh->mesh_map("MATERIAL");

  for (int i = 0; i < d_precursors_group; ++i)
  {
    double inv_lambda = 1.0 / d_material->lambda(i);
    for (int cell = 0; cell < d_num_cells; ++cell)
   {
     (*d_P0)[cell + i*d_num_cells] = inv_lambda * d_material->beta(mt[cell], i) * fd[cell];
     // store the initial concentration in the solution matrix
     (*d_precursors)(cell + i*d_num_cells, 0) =  (*d_P0)[cell + i*d_num_cells];
   }
 }
}

//------------------------------------------------------------------------------------//

/// Project the initial flux and precursors concentraion onto their basis to get
void TransientSolver::ProjectInitial()
{
  for (int g=0; g<d_number_groups; g++)
  {
    for (int cell=0; cell<d_num_cells; cell++)
    {
      (*d_phi0)[cell + g*d_num_cells] = d_state->phi(g)[cell];
      //store the initial flux
      (*d_flux)(cell + g*d_num_cells, 0) = d_state->phi(g)[cell];
    }
  }

  d_flux_basis->multiply_transpose(*d_phi0, *d_phi0_r);

  // project initial precursors
  d_precursors_basis->multiply_transpose(*d_P0, *d_P0_r);

  // stack the flux and precursors in one vector
  for (int i =0; i<d_rf; i++)
  {
    (*d_sol0_r)[i] = (*d_phi0_r)[i];
  }

  for (int i =0; i<d_rc; i++)
  {
    (*d_sol0_r)[d_rf + i] = (*d_P0_r)[i];
  }
}

//------------------------------------------------------------------------------------//

void TransientSolver::Construct_Operator(double t, double dt)
{
  // get the matrices
  KineticMatrices K(d_inp, d_mesh, d_material, d_flux_basis, d_precursors_basis);

  d_precursors_decay = new callow::MatrixDense(d_rc, d_rc);
  d_precursors_decay = K.precursors_decay();

  d_delayed_production = new callow::MatrixDense(d_rf, d_rc);
  d_delayed_production = K.delayed_production();

  d_precursors_production = new callow::MatrixDense(d_rc, d_rf);
  d_precursors_production = K.precursors_production();

  // if diffusion
  SP_lossoperator L(new DiffusionLossOperator(d_inp, d_material, d_mesh, false, 0.0, false, 1.0, false));
  d_L = L;

  SP_gainoperator G(new DiffusionGainOperator(d_inp, d_material, d_mesh, false));
  d_G = G;
  //d_G = new DiffusionGainOperator(d_inp, d_material, d_mesh, false);

 callow::Matrix::SP_matrix d_L_prime;
 callow::Matrix::SP_matrix d_G_prime;

  int m =d_number_groups*d_num_cells;

  vec_int nnz_L(m, 1 + 2 * d_mesh->dimension() + d_number_groups);
  d_L_prime = new callow::Matrix(d_number_groups*d_num_cells, d_number_groups*d_num_cells, nnz_L[0]);

  vec_int nnz_G(m, d_number_groups);
  d_G_prime = new callow::Matrix(d_number_groups*d_num_cells, d_number_groups*d_num_cells, nnz_G[0]);

  int * rows_L = d_L->rows();
  int* cols_L = d_L->columns();
  double* v_L = d_L->values();

  int * rows_G = d_G->rows();
  int* cols_G = d_G->columns();
  double* v_G = d_G->values();

  OperatorProjection Projector(1);
  d_Gr = new callow::MatrixDense(d_rf, d_rf);
  d_Lr = new callow::MatrixDense(d_rf, d_rf);


  const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");

  for (int g=0; g < d_number_groups; g++)
  {
    // loop over cells
    for (int cell=0; cell < d_num_cells; cell++)
    {
      int m = mat_map[cell];
      int r = d_num_cells*g + cell; //row number
      for (int p = rows_L[r]; p < rows_L[r + 1]; p++)
      {
	int c = cols_L[p];
        double value = -v_L[p]*d_material->velocity(g);
        //d_L->insert(r, c, value, 0);
        d_L_prime->insert(r, c, value);
      }

      for (int p = rows_G[r]; p < rows_G[r + 1]; p++)
      {
        int c = cols_G[p];
        double value = v_G[p]*(1 - d_material->beta_total(m))*d_material->velocity(g);
        //d_G->insert(r, c, value, 0);
        d_G_prime->insert(r, c, value);
      }
    }
 }

 d_L_prime->assemble();
 d_G_prime->assemble();

  Projector.SetOperators(d_G_prime, d_flux_basis);

  Projector.Project(d_Gr);

  Projector.SetOperators(d_L_prime, d_flux_basis);
  Projector.Project(d_Lr);



  // this assumes that a basis set is generated for each flux group, so the velocity
  // is not collapsed.
  int r = d_rf/d_number_groups;
  for (int g=0; g<d_number_groups; g++)
  {
    for (int i=0; i<r; i++)
    { // fill the upper left of A
      for (int j=0; j<d_rf; j++)
      {
        (*d_A)(i + g*r, j) = ((*d_Lr)(i + g*r, j) + (*d_Gr)(i + g*r, j));
      }
      // fill the upper right
      for (int k =0; k<d_rc; k++)
      {
        (*d_A)(i + g*r, k+d_rf) = (*d_delayed_production)(i+g*r, k)*d_material->velocity(g);
      }
    }
  }
   // fill the lower half
  for (int i=0; i<d_rc ; i++)
  {  // lower left
     for (int j=0; j<d_rf; j++)
     {
        (*d_A)(d_rf+i, j) = (*d_precursors_production)(i, j);
     }
    // lower right
    for (int k=0; k<d_rc; k++)
    {
      (*d_A)(d_rf+i, k+d_rf) = (*d_precursors_decay)(i, k);
    }
  }
 }

//---------------------------------------------------------------------//

void TransientSolver::Refersh_Operator()
{
  // if diffusion
  SP_lossoperator L(new DiffusionLossOperator(d_inp, d_material, d_mesh, false, 0.0, false, 1.0, false));
  d_L = L;

  int * rows_L = d_L->rows();
  int* cols_L = d_L->columns();
  double* v_L = d_L->values();
  vectorized_matrix = v_L;

  if (deim_flag)
  {
   online();
  }
  else
  {
  d_Lr = new callow::MatrixDense(d_rf, d_rf);

  OperatorProjection Projector(1);
  Projector.SetOperators(d_L, d_flux_basis);
  Projector.Project(d_Lr);
  }
  const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");

  // this assumes that a basis set is generated for each flux group, so the velocity
  // is not collapsed.
  int r = d_rf/d_number_groups;

  for (int g=0; g<d_number_groups; g++)
  {
    for (int i=0; i<r; i++)
    {
      for (int j=0; j<d_rf; j++)
      {
        (*d_A)(i + g*r, j) = (-(*d_Lr)(i + g*r, j)*d_material->velocity(g) + (*d_Gr)(i + g*r, j));
      }
    }
  }
}
//------------------------------------------------------------------------------------//

void TransientSolver::Solve(SP_state initial_state)
{

  d_material->update(0.0, 0, 1, false);

  Construct_Operator(0, d_dt);

  d_state = initial_state;

  initialize_precursors();

  ProjectInitial();

 // if (d_multiphysics) *d_vec_multiphysics[0] = *d_multiphysics;

  double t = 0.0;

  int n = d_rf + d_rc;

  SP_matrix d_A_;

  d_A_ = new callow::MatrixDense(n, n);

  if (deim_flag)
  {
    offline();
  }

  for (int step=0 ; step< d_number_steps; step++)
  {
	std::cout << "step  = " << step << "\n";

    t += d_dt;

    d_material->update(t, d_dt, 1, false);

    Refersh_Operator();

    for (int i=0; i<n; i++)
    {
      for (int j=0; j<n; j++)
      {
        (*d_A_)(i, j) = -(*d_A)(i, j)*d_dt;
        if (i == j) (*d_A_)(i, j) += 1;
      }
    }

   d_solver->set_operators(d_A_, db);

   d_solver->solve(*d_sol0_r, *d_sol_r);


   if (step == 0)
   {
   	 d_sol0_r->print_matlab("b.txt");
   	 d_sol_r->print_matlab("x.txt");
     d_A_->print_matlab("A.txt");
   }

   // store this state in the solution matrix
   for (int i=0; i<n; i++)
   {
     //(*d_sols)(i, step) = (*d_sol_r)[i];
     (*d_sol0_r)[i] = (*d_sol_r)[i];
   }

   // reconstruct the full solution
   reconstruct(step);

   if (d_multiphysics)
   {
     update_multiphysics(t, d_dt, 1);
     //*d_multiphysics_0 = *d_multiphysics;
     *d_vec_multiphysics[0] = *d_multiphysics;
   }

  }
  // temporary ... need to have getter for flux, etc
  d_flux->print_matlab("flux.txt");
  d_precursors->print_matlab("precursors.txt");
}

//------------------------------------------------------------------------------------//

void TransientSolver::reconstruct(int step)

{
  callow::Vector phi(d_number_groups*d_num_cells, 0.0);
  callow::Vector C(d_precursors_group*d_num_cells, 0.0);

  callow::Vector v1(d_rf, 0.0);
  callow::Vector v2(d_rc, 0.0);

  for (int f=0; f< d_rf; f++)
  {
    v1[f] = (*d_sol_r)(f);
  }
  d_flux_basis->multiply(v1, phi);

  for (int p=0; p<d_rc; p++)
  {
    v2[p] = (*d_sol_r)(p + d_rf);
  }

  d_precursors_basis->multiply(v2, C);
  double *phi_ = &phi[0];
  double *C_ = &C[0];

  d_flux->insert_col(step+1, phi_, 0);
  d_precursors->insert_col(step+1, C_, 0);

  //std::cout << "phi[0][0]= " << phi[0] << " \n";

  ///
  fluxes.clear();

  for (int g=0; g<d_number_groups; g++)
  {
    callow::Vector phi_g(d_num_cells, 0.0);
    for (int i=0; i<d_num_cells; i++)
    {
     phi_g[i] = phi[d_num_cells*g + i];
    }

    fluxes.push_back(phi_g);
  }
 // std::cout << "fluxes[0][0]= " << fluxes[0][0] << " ***\n";
  //std::cout << "reconstructed &&&&&&&&&&&&&&&&&&\n";

}

//------------------------------------------------------------------------------//

void TransientSolver::offline()
{
   DEIM D(d_deim_basis, 3);
   D.Search();
   l = D.interpolation_indices();
   Ur_deim = new callow::MatrixDense(3, 3);
   Ur_deim = D.ReducedBasis();

   offline_stage O(d_L, d_flux_basis, d_deim_basis, 3);

   M_L = O.Decompositon();

   LinearSolverCreator::SP_db db_deim;
   db_deim = get_db();

   d_solver_deim = LinearSolverCreator::Create(db_deim);

   d_solver_deim->set_operators(Ur_deim, db_deim);
}

//----------------------------------------------------------------------------//
void TransientSolver::online()
{

	d_x_deim = new callow::Vector(3, 0.0);
	d_b_deim = new callow::Vector(3, 0.0);

	for (int i=0; i<3; i++)
	{
	  int id = l[i];
     (*d_b_deim)(i) = vectorized_matrix[id];
	}

	d_solver_deim->solve(*d_b_deim, *d_x_deim);

	d_Lr = new callow::MatrixDense(d_rf, d_rf);

	double v;
	// construct the Left operator
	for (int r=0; r< 3; r++)
	{
	  for (int i=0; i<d_rf; i++)
	  {
	    for (int j=0; j<d_rf; j++)

	  {
	    v = (*d_x_deim)[r]*(*M_L[r])(i, j);
	    d_Lr->insert(i, j, v, 1);
	  }
	 }
	}

}

//----------------------------------------------------------------------------//

void TransientSolver::
set_multiphysics(SP_multiphysics ic,
                 multiphysics_pointer_rom update_multiphysics_rhs_rom,
                 void* multiphysics_data)
{

  int d_order= 1;
  // Preconditions
  Require(ic);
  Require(update_multiphysics_rhs_rom);

  std::cout << " ic T=" << ic->variable(0)[0] << std::endl;

  d_multiphysics            = ic;
  d_update_multiphysics_rhs_rom = update_multiphysics_rhs_rom;
  d_multiphysics_data       = multiphysics_data;

  // Create previous physics states
  d_vec_multiphysics.resize(d_order);
  d_multiphysics_0 = new MultiPhysics(*ic);
  for (int i = 0; i < d_order; ++i)
  {
    d_vec_multiphysics[i] = new MultiPhysics(*ic);
  }
  std::cout << d_vec_multiphysics[0]->variable(0)[0] << std::endl;
  std::cout << " ic T=" << d_multiphysics->variable(0)[0] << std::endl;
}

//----------------------------------------------------------------------------//

void TransientSolver::
update_multiphysics(const double t, const double dt, const size_t order)
{
  std::cout << "t_eval " << t << "\n";
  // Update the right hand side.  The result is placed into
  // the working vector d_multiphysics
  //std::cout << " P before = " << d_multiphysics->variable(0)[0] - 300.0 << std::endl;
  d_update_multiphysics_rhs_rom(d_multiphysics_data, fluxes, t, dt);
  //std::cout << " P after = " << d_multiphysics->variable(0)[0] << std::endl;

  // Loop through and compute
  //  y(n+1) = (1/a0) * ( dt*rhs + sum of bdf terms )
  for (size_t i = 0; i < d_multiphysics->number_variables(); ++i)
  {
    // Reference to P(n+1)
    MultiPhysics::vec_dbl &P   = d_multiphysics->variable(i);

    //std::cout << " Pold[0]=" << P[0] << std::endl;
    //printf("delP = %18.12e \n", P[0]);
   // printf("Pold[0] = %18.12e \n", d_vec_multiphysics[0]->variable(0)[0]);

    // Loop over all elements (usually spatial)
    for (int j = 0; j < P.size(); ++j)
    {
      double v = dt * P[j];
      for (size_t k = 1; k <= order; ++k)
        v += bdf_coefs[order-1][k] * d_vec_multiphysics[k-1]->variable(i)[j];
      P[j] = v / bdf_coefs[order-1][0];
    } // end element loop
    std::cout << " Pnew[0]=" << P[0] - 300.0 << std::endl;
  } // end variable loop
}

