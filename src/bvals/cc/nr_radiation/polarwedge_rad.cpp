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

// Helper: given index of angle n, returns index of angle achieved by
//         flipping angle in both theta/phi -- necessary for corotating frame
int FlipAngleAcrossPoleCorotating(MeshBlock *pmb, int n) {
  // assumes 2D spherical polar

  NRRadiation *prad = pmb->pnrrad;
  const int nang = prad->nang;
  int flipped_n;

  // ANGLE_FLAG_COROTATING makes a "corotating" angle grid whose basis vectors
  // are the basis vectors of the coordinate system
  //
  // Accomplish a flip by flipping the angles according to what is appropriate
  //  angles are indexed as:
  //    n = 2*npsi*coszeta_idx + psi_idx
  //
  // Angles have LOCAL components:
  //  x = sinzeta cos(psi)
  //  y = sinzeta sin(psi)
  //  z = coszeta
  const int npsi = prad->npsi;

  // special case of polar_angle -- polar rays left unchanged
  if (prad->polar_angle && n >= nang-2) return n;

  int coszeta_idx = n / (2*npsi);
  int psi_idx = n - 2*npsi*coszeta_idx;
  
  int flipped_psi_idx = (psi_idx + npsi) % (2*npsi);
  flipped_n = 2*npsi*coszeta_idx + flipped_psi_idx;

  return flipped_n;
}

// pmb: MeshBlock being updated; var: its specific intensities (*var_cc)
// face_idx: jl (inner) or ju (outer); is_outer: 0=inner 1=outer
// copy specific intensity into ghost cells, applying flip in theta/phi if necessary
void ApplyPolarWedge(MeshBlock *pmb, AthenaArray<Real> &var,
                     int il, int iu, int face_idx, int kl, int ku, int ngh,
                     int is_outer) {
  const int& nang = pmb->pnrrad->nang; // number of angles
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency groups
  const int& n_fre_ang = pmb->pnrrad->n_fre_ang; // nang*nfreq

  for (int k=kl; k<=ku; ++k) { // loop over spatial grid (ghost cells in theta -> j)
    for (int j=1; j<=ngh; ++j) {

      // depending on inner/outer, pick ghost/active cell as j'th cell from face_idx
      int jg, ja;
      if (!is_outer) { jg = face_idx-j;  ja = face_idx+j-1; }
      if ( is_outer) { jg = face_idx+j;  ja = face_idx-j+1; }

      for (int i=il; i<=iu; ++i) {

        if (pmb->pnrrad->angle_flag == ANGLE_FLAG_COROTATING) {

          for (int n=0; n < nang; ++n) {
            // for corotating-frame angular grids, the theta/phi directions both swap
            int m = FlipAngleAcrossPoleCorotating(pmb, n);

            for (int ifr=0; ifr<nfreq; ++ifr) { // loop over frequency groups
              Real *ir_active = &var(k, ja, i, ifr*nang + m);  // active cell, flipped angle m
              Real *ir_ghost  = &var(k, jg, i, ifr*nang + n);  // ghost cell, angle n
              *ir_ghost = *ir_active;
            }
          }
        } else {
          // for lab-frame angular grids, just copy ENTIRE specific intensities across
          // loop over all frequency groups / angles
          #pragma omp simd // attempt optimization
          for (int n=0; n<n_fre_ang; ++n) var(k,jg,i,n) = var(k,ja,i,n);
        }

      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::PolarWedgeInnerX2(
//          Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//  \brief polar wedge boundary conditions, inner x2 boundary
void RadBoundaryVariable::PolarWedgeInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  ApplyPolarWedge(pmy_block_, *var_cc, il, iu, jl, kl, ku, ngh, 0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::PolarWedgeOuterX2(
//          Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//  \brief polar wedge boundary conditions, outer x2 boundary
void RadBoundaryVariable::PolarWedgeOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  ApplyPolarWedge(pmy_block_, *var_cc, il, iu, ju, kl, ku, ngh, 1);
}
