//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundaryTraits.hh
 *  @author robertsj
 *  @date   Sep 11, 2012
 *  @brief  BoundaryTraits class definition.
 *  @todo   Specialization for MOC needed
 */
//---------------------------------------------------------------------------//

#ifndef detran_BOUNDARYTRAITS_HH_
#define detran_BOUNDARYTRAITS_HH_

#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @brief Boundary traits to simplify type access.
 *
 *  For fine mesh discretizations, the boundary data is stored for
 *  each surface in the form of a value (1-D), 1-D vector (2-D), or
 *  2-D vector (3-D).  This would be used for mesh-based discretizations,
 *  as in SN and diffusion.
 */
//---------------------------------------------------------------------------//

template <class D>
struct BoundaryTraits
{
  /* ... */
};

template <>
struct BoundaryTraits<_3D>
{
  typedef detran_utilities::vec2_dbl value_type;
};

template <>
struct BoundaryTraits<_2D>
{
  typedef detran_utilities::vec_dbl value_type;
};

template <>
struct BoundaryTraits<_1D>
{
  typedef double value_type;
};

//---------------------------------------------------------------------------//
/**
 *  @brief Boundary access.
 *
 *  Because we employ a templated boundary type, it can sometimes
 *  complicate simple access within an otherwise general algorithm.
 *  This accessor returns a boundary element given boundary spatial
 *  indices.
 */
//---------------------------------------------------------------------------//

template <class D>
struct BoundaryValue
{
  /* ... */
};

template <>
struct BoundaryValue<_3D>
{
  // Mutable access to boundary value
  static inline double&
  value(BoundaryTraits<_3D>::value_type &b,
        const size_t                     i,
        const size_t                     j)
  {
    Require(i < b.size());
    Require(j < b[0].size());
    return b[i][j];
  }
};

template <>
struct BoundaryValue<_2D>
{
  // Mutable access to boundary value
  static inline double&
  value(BoundaryTraits<_2D>::value_type &b,
        const size_t                     i,
        const size_t                     j = 0)
  {
    Require(i < b.size());
    return b[i];
  }
};

template <>
struct BoundaryValue<_1D>
{
  // Mutable access to boundary value
  static inline double&
  value(BoundaryTraits<_1D>::value_type &b,
        const size_t                     i = 0,
        const size_t                     j = 0)
  {
    return b;
  }
};

} // end namespace detran

#endif /* detran_BOUNDARYTRAITS_HH_ */
