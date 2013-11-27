//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_DD_3D.i.hh
 *  @author Jeremy Roberts
 *  @date   Mar 31, 2012
 *  @brief  Equation_DD_3D inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_DD_3D_I_HH_
#define detran_EQUATION_DD_3D_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
inline void Equation_DD_3D::solve(const size_t i,
                                  const size_t j,
                                  const size_t k,
                                  moments_type &source,
                                  face_flux_type &psi_in,
                                  face_flux_type &psi_out,
                                  moments_type &phi,
                                  angular_flux_type &psi)
{
  typedef detran_geometry::CartesianMesh Mesh;

  // Compute cell-center angular flux.
  int cell = d_mesh->index(i, j, k);
  double coef = 1.0 / (d_material->sigma_t(d_mat_map[cell], d_g) +
                       d_coef_x[i] + d_coef_y[j] + d_coef_z[k]);
  double psi_center = coef * (source[cell] + d_coef_x[i] * psi_in[Mesh::YZ] +
                                             d_coef_y[j] * psi_in[Mesh::XZ] +
                                             d_coef_z[k] * psi_in[Mesh::XY]);

  // Compute outgoing fluxes.
  double two_psi_center = 2.0 * psi_center;
  psi_out[Mesh::YZ] = two_psi_center - psi_in[Mesh::YZ];
  psi_out[Mesh::XZ] = two_psi_center - psi_in[Mesh::XZ];
  psi_out[Mesh::XY] = two_psi_center - psi_in[Mesh::XY];

  // Compute flux moments.
  phi[cell] += d_quadrature->weight(d_angle) * psi_center;

  // Store angular flux if needed.
  if (d_update_psi) psi[cell] = psi_center;

}

} // end namespace detran

#endif /* detran_EQUATION_DD_3D_I_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_DD_3D.i.hh
//---------------------------------------------------------------------------//
