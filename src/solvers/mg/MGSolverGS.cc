//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGSolverGS.cc
 *  @author robertsj
 *  @date   Apr 10, 2012
 *  @brief  MGSolverGS member definitions.
 */
//---------------------------------------------------------------------------//

#include "MGSolverGS.hh"
#include <iostream>
#include <string>

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
MGSolverGS<D>::MGSolverGS(SP_state                  state,
                          SP_material               material,
                          SP_boundary               boundary,
                          const vec_externalsource &q_e,
                          SP_fissionsource          q_f,
                          bool                      multiply)
  : Base(state, material, boundary, q_e, q_f, multiply)
  , d_norm_type("Linf")
{
  if (d_input->check("outer_norm_type"))
    d_norm_type = d_input->template get<std::string>("outer_norm_type");

  // Post conditions
  Ensure(d_norm_type == "Linf" || d_norm_type == "L1" || d_norm_type == "L2");
}

//---------------------------------------------------------------------------//
template <class D>
int MGSolverGS<D>::number_sweeps() const
{
  return d_wg_solver->get_sweeper()->number_sweeps();
}

template class MGSolverGS<_1D>;
template class MGSolverGS<_2D>;
template class MGSolverGS<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of GaussSeidelMG.cc
//---------------------------------------------------------------------------//

