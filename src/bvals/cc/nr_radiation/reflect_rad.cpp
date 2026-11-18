//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reflect.cpp
//  \brief implementation of reflecting BCs in each dimension

// C headers

// C++ headers
// #include <cstring>  // std::strcmp, for the commented-out polar_angle fix

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../coordinates/coordinates.hpp"
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

int MirrorAngle(MeshBlock *pmb, int n, int dim_face);

// pmb: MeshBlock being updated; var: its specific intensities (*var_cc)
void Reflect(MeshBlock *pmb, AthenaArray<Real> &var,
             int il, int iu, int jl, int ju, int kl, int ku, int ngh,
             int is_outer, int dim_face) {
  const int& nang = pmb->pnrrad->nang; // number of angles
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency groups

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

          for (int n=0; n < nang; ++n) { // loop over angles
            int m = MirrorAngle(pmb, n, dim_face);  // reflected angle index, depends on angle_flag
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

// Helper: given index of angle n, returns index of "mirror" angle corresponding
//         to flipping angle n across dim_face ("reflection")
int MirrorAngle(MeshBlock *pmb, int n, int dim_face) {
  NRRadiation *prad = pmb->pnrrad;
  const int nang = prad->nang;
  const int angle_flag = prad->angle_flag;
  int mirror_n;

  if (angle_flag==ANGLE_FLAG_LAB_DAVIS || angle_flag==ANGLE_FLAG_LAB_MAGIC) {
    // WARNING (NZR): angle_flag=0,10 may not function correctly in curvilinear coordinates
    //   this is because octant flipping is not the correct way to impose reflections in
    //   coordinates

    // angles in the same "octant" are adjacent -- identify octant
    // our angle is n = oct * nang_per_oct + ang_in_oct
    int nang_per_oct = nang / prad->noct;
    int oct = n / nang_per_oct; // which octant
    int ang_in_oct = n % nang_per_oct; // which angle in octant

    // octants are indexed like this:
    //  1D: oct = 1*pos(mu_x)
    //  2D: oct = 2*pos(mu_y) + 1*pos(mu_x)
    //  3D: oct = 4*pos(mu_z) + 2*pos(mu_y) + 1*pos(mu_x)
    // where pos(mu) = 0(mu positive) and 1(mu negative)
    
    // for 2D/3D, get mirror octant by flipping bit of face dimension
    int mirror_oct = oct ^ (1 << dim_face);
    mirror_n = mirror_oct * nang_per_oct + ang_in_oct;

  } else if (angle_flag == ANGLE_FLAG_COROTATING) {
    // ANGLE_FLAG_COROTATING makes a "corotating" angle grid whose basis vectors
    // are the basis vectors of the coordinate system
    //
    // Accomplish a mirror by flipping the angles according to what is appropriate
    //  angles are indexed as:
    //    n = 2*npsi*coszeta_idx + psi_idx
    //
    // Angles have LOCAL components:
    //  x = sinzeta cos(psi)
    //  y = sinzeta sin(psi)
    //  z = coszeta
    // But at some face (say x) the local x grid coincides with the global x one,
    // so flipping the angular index in x produces a reflection in x.
    const int nzeta = prad->nzeta;
    const int npsi = prad->npsi;

    // coordinate direction (0, 1, 2 = x1, x2, x3) of the local x, y, z axes;
    // spherical_polar: (axisx, axisy, axisz) = (1, 2, 0)
    int axisx = 0, axisy = 1, axisz = 2;
    pmb->pcoord->AxisDirection(&axisx, &axisy, &axisz);

    // special case of polar_angle -- last two angles are pole rays
    //   if reflection is z, swap the poles. otherwise, leave it the same
    //   (poles are invariant under x/y reflections)
    if (prad->polar_angle && n >= nang-2) {
      // locate dimension corresponding to added polar angles
      int ndim = 1;
      if (pmb->block_size.nx2 > 1) ndim = 2;
      if (pmb->block_size.nx3 > 1) ndim = 3;
      bool spherical_polar = (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0);

      int dim_pole;
      if (ndim == 1) {
        dim_pole = 0;
      } else if (ndim == 2 && !spherical_polar && npsi > 1) {
        // poles are along axisx only for the 2D non-spherical-polar psi grid;
        // with npsi == 1 they are along axisz (the default below)
        dim_pole = axisx;
      } else {
        dim_pole = axisz;
      }

      if (dim_pole == dim_face) { // if reflection along pole, switch poles
        if (n == nang-2) {
          mirror_n = nang-1;
        } else {
          mirror_n = nang-2;
        }
      } else { // if reflection not along pole, pole left unchanged
        mirror_n = n;
      }

      return mirror_n;
    }

    // special case of 1D -- there is no psi grid (npsi = 0), angles are indexed
    //   n = coszeta_idx, and all rays lie along x1
    if (pmb->block_size.nx2 == 1) {
      if (dim_face == 0) { // flip coszeta
        mirror_n = 2*nzeta - 1 - n;
      } else {
        mirror_n = n;
      }
      return mirror_n;
    }

    int coszeta_idx = n / (2*npsi);
    int psi_idx = n - 2*npsi*coszeta_idx;
    int mirror_coszeta_idx, mirror_psi_idx;

    // special case of 2D with npsi == 1 -- the two psi indices are at psi = 0, pi
    //   (not pi/2, 3pi/2 as the general psi grid would give), so the rays have
    //   x = +/- sinzeta and y = 0
    bool is_2d_npsi1 = (pmb->block_size.nx3 == 1 && npsi == 1);

    if (dim_face == axisx) {
      // flip cos(psi) while preserving sin(psi): psi -> pi - psi, mod 2pi
      //  add another 2pi because % operator in C++ does not implement true
      //  modulo for negative numbers
      mirror_coszeta_idx = coszeta_idx;
      if (is_2d_npsi1) { // swap psi = 0 and psi = pi
        mirror_psi_idx = 1 - psi_idx;
      } else {
        mirror_psi_idx = (3*npsi - 1 - psi_idx) % (2*npsi);
      }

    } else if (dim_face == axisy) {
      // flip sin(psi) while preserving cos(psi): psi -> 2pi - psi
      mirror_coszeta_idx = coszeta_idx;
      if (is_2d_npsi1) { // y = 0, so the ray is unchanged
        mirror_psi_idx = psi_idx;
      } else {
        mirror_psi_idx = 2*npsi - 1 - psi_idx;
      }

    } else if (dim_face == axisz) {
      // flip coszeta
      mirror_coszeta_idx = 2*nzeta - 1 - coszeta_idx;
      mirror_psi_idx = psi_idx;

    }

    mirror_n = 2*npsi*mirror_coszeta_idx + mirror_psi_idx;
  }

  return mirror_n;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x1 boundary
void RadBoundaryVariable::ReflectInnerX1(Real time, Real dt, int il, int jl, int ju,
                                         int kl, int ku, int ngh) {
  Reflect(pmy_block_, *var_cc, il,  0, jl, ju, kl, ku, ngh, 0, 0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX1(
//          Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x1 boundary
void RadBoundaryVariable::ReflectOuterX1(Real time, Real dt, int iu, int jl, int ju,
                                         int kl, int ku, int ngh) {
  Reflect(pmy_block_, *var_cc,  0, iu, jl, ju, kl, ku, ngh, 1, 0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX2(
//          Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x2 boundary
void RadBoundaryVariable::ReflectInnerX2(Real time, Real dt, int il, int iu, int jl,
                                         int kl, int ku, int ngh) {
  Reflect(pmy_block_, *var_cc, il, iu, jl,  0, kl, ku, ngh, 0, 1);
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX2(
//          Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x2 boundary
void RadBoundaryVariable::ReflectOuterX2(Real time, Real dt, int il, int iu, int ju,
                                         int kl, int ku, int ngh) {
  Reflect(pmy_block_, *var_cc, il, iu,  0, ju, kl, ku, ngh, 1, 1);
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh)
//  \brief REFLECTING boundary conditions, inner x3 boundary
void RadBoundaryVariable::ReflectInnerX3(Real time, Real dt, int il, int iu, int jl,
                                         int ju, int kl, int ngh) {
  Reflect(pmy_block_, *var_cc, il, iu, jl, ju, kl,  0, ngh, 0, 2);
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x3 boundary
void RadBoundaryVariable::ReflectOuterX3(Real time, Real dt, int il, int iu, int jl,
                                         int ju, int ku, int ngh) {
  Reflect(pmy_block_, *var_cc, il, iu, jl, ju,  0, ku, ngh, 1, 2);
}
