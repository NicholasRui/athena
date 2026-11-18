//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reflect.cpp
//  \brief implementation of reflecting BCs in each dimension

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../nr_radiation/radiation.hpp"
#include "./bvals_rad.hpp"

// At each face (e.g., defined by constant X), loop over the entire "face" of
//   ghost cells, i.e., loop over "thickness" of ghost cell layer in X and the full area
//   Y and Z. Each ghost cell lies some distance S within the "ghost region."
//  
// To implement reflection, we set the specific intensities of each ghost cell (lying
//   a distance S within the ghost region) by copying the specific intensities of a
//   "partner" active cell (lying a distance S within the active region) but with
//   the angle grid flipped along that axis.
//
// PSEUDOCODE: to reflect across X...
  // loop over ghost cells X_ghost=1...S, Y[active], Z[active]:
  //   loop over ghost angles n_ghost:
  //     n_active = flip n_ghost across X
  //     X_active = X[active cell a distance X_ghost away from face]
  //     I[n_ghost, (X, Y, Z)] = I[n_active, (X_active, Y, Z)]
//
// There are six different possibilities, depending on dimension X/Y/Z and which end of
//   the domain inner/outer
//
// Define a helper that executes this:
void RadBoundaryVariable::Reflect(int il, int iu, int jl, int ju, int kl, int ku, int ngh,
                                  int is_outer, int dim_face) {
  const int& nang = pmy_block_->pnrrad->nang; // number of angles
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency groups
  AthenaArray<Real> &var = *var_cc;

  // is_outer: 0=inner 1=outer
  // dim_face: 0=x, 1=y, 2=z
  int face_idx, pl, pu, ql, qu;
  if (dim_face == 0 && !is_outer) { face_idx = il; pl = jl; pu = ju; ql = kl; qu = ku; }
  if (dim_face == 0 &&  is_outer) { face_idx = iu; pl = jl; pu = ju; ql = kl; qu = ku; }
  if (dim_face == 1 && !is_outer) { face_idx = jl; pl = il; pu = iu; ql = kl; qu = ku; }
  if (dim_face == 1 &&  is_outer) { face_idx = ju; pl = il; pu = iu; ql = kl; qu = ku; }
  if (dim_face == 2 && !is_outer) { face_idx = kl; pl = il; pu = iu; ql = jl; qu = ju; }
  if (dim_face == 2 &&  is_outer) { face_idx = ku; pl = il; pu = iu; ql = jl; qu = ju; }

  for (int ifr=0; ifr<nfreq; ++ifr) { // loop over frequency groups

    for (int p=pl; p<=pu; ++p) { // loop over perpendicular directions
      for (int q=ql; q<=qu; ++q) {
        for (int r=1; r<=ngh; ++r) { // loop over ghost direction

          // depending on inner/outer, pick ghost/active cell as r'th cell from face_idx
          int ghost_face_idx, active_face_idx;
          if (!is_outer) { ghost_face_idx = face_idx-r;  active_face_idx = face_idx+r-1; }
          if ( is_outer) { ghost_face_idx = face_idx+r;  active_face_idx = face_idx-r+1; }

          // based on dim_face, find indices for partner ghost (ig,jg,kg) and active (ia,ja,ka) cells
          int ig, jg, kg, ia, ja, ka;
          if (dim_face == 0) { ig=ghost_face_idx; ia=active_face_idx; jg=ja=p; kg=ka=q; }
          if (dim_face == 1) { jg=ghost_face_idx; ja=active_face_idx; ig=ia=p; kg=ka=q; }
          if (dim_face == 2) { kg=ghost_face_idx; ka=active_face_idx; ig=ia=p; jg=ja=q; }

          for (int n = 0; n < nang; ++n) { // loop over angles
            int m = MirrorAngle(n, dim_face);  // reflected angle index, depends on angle_flag
            Real *ir_active = &var(ka, ja, ia, ifr*nang + m);  // active cell, reflected angle m
            Real *ir_ghost  = &var(kg, jg, ig, ifr*nang + n);  // ghost cell, angle n
            *ir_ghost = *ir_active;
          }

        }
      }
    }
  }

  return;
}


// TODO: need MirrorAngle


void RadBoundaryVariable::MirrorAngle(int n, int dim_face) {
  int angle_flag = pmy_block_->pnrrad->angle_flag;

  if (angle_flag==ANGLE_FLAG_LAB_DAVIS || angle_flag==ANGLE_FLAG_LAB_MAGIC) {



  // TODO!!








  } else if (angle_flag == ANGLE_FLAG_COROTATING) {

  }
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x1 boundary
void RadBoundaryVariable::ReflectInnerX1(Real time, Real dt, int il, int jl, int ju,
                                         int kl, int ku, int ngh) {
  Reflect(il,  0, jl, ju, kl, ku, ngh, 0, 0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX1(
//          Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x1 boundary
void RadBoundaryVariable::ReflectOuterX1(Real time, Real dt, int iu, int jl, int ju,
                                         int kl, int ku, int ngh) {
  Reflect( 0, iu, jl, ju, kl, ku, ngh, 1, 0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX2(
//          Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x2 boundary
void RadBoundaryVariable::ReflectInnerX2(Real time, Real dt, int il, int iu, int jl,
                                         int kl, int ku, int ngh) {
  Reflect(il, iu, jl,  0, kl, ku, ngh, 0, 1);
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX2(
//          Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x2 boundary
void RadBoundaryVariable::ReflectOuterX2(Real time, Real dt, int il, int iu, int ju,
                                         int kl, int ku, int ngh) {
  Reflect(il, iu,  0, ju, kl, ku, ngh, 1, 1);
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh)
//  \brief REFLECTING boundary conditions, inner x3 boundary
void RadBoundaryVariable::ReflectInnerX3(Real time, Real dt, int il, int iu, int jl,
                                         int ju, int kl, int ngh) {
  Reflect(il, iu, jl, ju, kl,  0, ngh, 0, 2);
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x3 boundary
void RadBoundaryVariable::ReflectOuterX3(Real time, Real dt, int il, int iu, int jl,
                                         int ju, int ku, int ngh) {
  Reflect(il, iu, jl, ju,  0, ku, ngh, 1, 2);
}
