//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SweepSource.i.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  SweepSource inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPSOURCE_I_HH_
#define SWEEPSOURCE_I_HH_

#include <iostream>

namespace detran
{
  template <class D>
  inline void SweepSource<D>::build_fixed(int g)
  {

    // Zero out moments source.
    d_fixed_group_source.assign(d_fixed_group_source.size(), 0.0);

    // Add external sources, if present.
    for (int i = 0; i < d_moment_external_sources.size(); i++)
    {
      for (int cell = 0; cell < d_mesh->number_cells(); cell++)
      {
        d_fixed_group_source[cell] +=
            d_moment_external_sources[i]->source(cell, g);
      }
    }

    // Add fission source, if present \todo
    if (d_fissionsource)
    {
      for (int cell = 0; cell < d_mesh->number_cells(); cell++)
      {
        //d_fixed_group_source[cell] = d_fissionsource;
      }
    }

  }

  template <class D>
  inline void SweepSource<D>::build_fixed_with_scatter(int g)
  {

    // Add the external and/or fission source first.
    build_fixed(g);

    // Add the in-scatter.
    d_scattersource->build_in_scatter_source(g, d_fixed_group_source);

  }

  template <class D>
  inline void SweepSource<D>::build_within_group_scatter(int g, const moments_type &phi)
  {
    // Zero out moments source.
    d_scatter_group_source.assign(phi.size(), 0.0);

    // Build within-group scattering
    d_scattersource->build_within_group_source(g, phi, d_scatter_group_source);



  }

  template <class D>
  inline void SweepSource<D>::build_total_scatter(int g, const State::vec_moments_type &phi)
  {
    // Zero out moments source.
    d_scatter_group_source.assign(phi.size(), 0.0);

    // Build within-group scattering
    d_scattersource->build_total_group_source(g, phi, d_scatter_group_source);

  }

  template <class D>
  inline const typename SweepSource<D>::sweep_source_type&
  SweepSource<D>::source(int g, int o, int a)
  {
    // Zero out.
    d_source.assign(d_source.size(), 0.0);

    // Add moment contributions.
    for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      // Add fixed contribution \todo need to use mtod
      d_source[cell] += d_fixed_group_source[cell] * inv_four_pi;
                       // (*d_MtoD)(o, a, 0, 0);

      // Add scatter contribution
      d_source[cell] += d_scatter_group_source[cell] * inv_four_pi;
                   //     (*d_MtoD)(o, a, 0, 0);

    }

    // Add discrete contributions if present.
    for (int i = 0; i < d_discrete_external_sources.size(); i++)
    {
      for (int cell = 0; cell < d_mesh->number_cells(); cell++)
      {
        d_source[cell] +=
            d_discrete_external_sources[i]->source(cell, g);
      }
    }
    return d_source;
  }

  template class SweepSource<_1D>;
  template class SweepSource<_2D>;
  template class SweepSource<_3D>;

} // namespace detran

#endif /* SWEEPSOURCE_I_HH_ */

//---------------------------------------------------------------------------//
//              end of SweepSource.i.hh
//---------------------------------------------------------------------------//
