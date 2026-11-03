//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file polarwedge_rad.cpp
//  \brief implementation of polar wedge radiation BCs in x2 direction

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../nr_radiation/radiation.hpp"
#include "bvals_rad.hpp"

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::PolarWedgeInnerX2(
//          Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//  \brief polar wedge boundary conditions, inner x2 boundary

// Claude: Polar wedge boundary: the x2 boundary sits at (or next to) the pole of a
// spherical polar grid that covers only a wedge in phi, so the cells across the pole
// are not in the domain. Ghost cell jl-j copies its mirror cell jl+j-1 for every ray,
// with no change of angle index ('do not flip the sign').
//
// Claude: That matches angle_flag = 0, where each index is a fixed Cartesian direction.
// BUG (angle_flag = 1): each index is a direction relative to the local (r, theta, phi)
// axes. Across the pole e_theta and e_phi both reverse (e_r does not), which is why the
// hydro polar wedge flips v_theta and v_phi (flip_across_pole_hydro). The radiation
// equivalent is to take the mirror cell's ray with (mu_r, -mu_theta, -mu_phi): psi ->
// psi + pi, i.e. ipsi -> (ipsi + npsi) mod 2*npsi with izeta unchanged. For npsi == 1
// that is ipsi 0 <-> 1, which matches AngularGrid's mu_theta = +/-sin(zeta) layout. The
// polar angles (along +/-r) stay as they are. Copying the same index instead turns rays
// heading toward the pole in the ghost cell into rays heading away from it.
// do not flip the sign of specific intensity
void RadBoundaryVariable::PolarWedgeInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
#pragma omp simd
        // Claude: n runs over all nu_+1 = n_fre_ang entries (every angle of every
        // frequency group) of the cell.
        for (int n=0; n<=nu_; ++n) {
          (*var_cc)(k,jl-j,i,n) = (*var_cc)(k,jl+j-1,i,n);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::PolarWedgeOuterX2(
//          Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//  \brief polar wedge boundary conditions, outer x2 boundary

// Claude: Outer x2 version: ghost cell ju+j copies mirror cell ju-j+1. Same angle_flag
// = 1 BUG.
void RadBoundaryVariable::PolarWedgeOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
#pragma omp simd
        for (int n=0; n<=nu_; ++n) {
          (*var_cc)(k,ju+j,i,n) = (*var_cc)(k,ju-j+1,i,n);
        }
      }
    }
  }
  return;
}
