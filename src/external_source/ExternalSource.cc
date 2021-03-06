/*
 * ExternalSource.cc
 *
 *  Created on: Apr 24, 2013
 *      Author: robertsj
 */

#include "ExternalSource.hh"

namespace detran_external_source
{

//----------------------------------------------------------------------------//
ExternalSource::ExternalSource(size_t        number_groups,
                               SP_mesh       mesh,
                               SP_quadrature quadrature,
                               bool          discrete)
  : d_mesh(mesh)
  , d_quadrature(quadrature)
  , d_number_groups(number_groups)
  , d_number_angles(-1)
  , d_discrete(discrete)
{
  // Preconditions
  Require(number_groups > 0);
  Require(mesh);

  // Set the angle count. We leave this in the constructor, since
  // the quadrature need not be set.
  if (d_quadrature) d_number_angles = d_quadrature->number_angles();
}

} // end namespace detran_external_source
