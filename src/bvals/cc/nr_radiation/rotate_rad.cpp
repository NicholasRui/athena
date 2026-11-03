//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file rotate_rad.cpp
//  \brief implementation of rotating radiation BCs in each dimension
//======================================================================================

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../nr_radiation/radiation.hpp"
#include "./bvals_rad.hpp"

// Claude: Overview. These 'rotating' boundaries are for simulating only a wedge of the
// sphere and filling its ghost cells from the other end of the domain by symmetry:
// rotate_phi = 1 for a pi/2 wedge in phi, rotate_phi = 2 for a pi wedge in phi, and
// rotate_theta = 1 at the inner x2 boundary. They run after the periodic copy has
// filled the ghost cells, and relabel the rays so that each ray's direction is turned
// along with the cell. Callers: ApplyPhysicalBoundaries (bvals.cpp),
// ApplyRadPhysicalBoundaries (bvals_rad.cpp) and the coarse-level pass in
// bvals_refine.cpp.
//
// Claude: Like reflect_rad.cpp, they work by permuting whole octant blocks, so they
// assume the angle_flag = 0 layout: fixed Cartesian directions stored octant by octant.
//
// Claude: BUG (angle_flag = 1): with angle_flag = 1 the rays are defined relative to
// the local r, theta, phi axes. A rotation about the polar axis carries the local axes
// at phi onto the local axes at phi + delta, so for rotate_phi the periodic copy alone
// is already correct, and these octant permutations scramble it. They also use
// nang/noct, which means nothing in the angle_flag = 1 layout (see reflect_rad.cpp).
// For angle_flag = 1 rotate_phi should do nothing to the angles, and rotate_theta needs
// its own mapping in (izeta, ipsi).
// The angular octant ( in x-y plane) is
//   1  |  0       5  |  4
//   -------      ---------
//   3  |  2       7  |  6
// in X-Z plane, it is
//   5  |  4       7  |  6
//   -------      ---------
//   1  |  0       3  |  2


// in radiation class, nang_per_oct is angles per octant, noct is the number of octant


// Claude: Rotation by pi about the z axis for one cell and one frequency group, in the
// z > 0 half only: (mu_x, mu_y) -> (-mu_x, -mu_y), i.e. octants 0 <-> 3 and 1 <-> 2.
// Unused: nothing calls CopyIntensity2. The 'direction' argument is ignored.
// Temporary function to copy intensity
// swap 0 and 3, 1 and 2
// change sign of both Fx and Fy
// only do this for 2D spherical polar
void CopyIntensity2(Real *ir, int nang_per_oct, int direction) {
  // here ir is only intensity for each cell and each frequency band
  for (int n=0; n<nang_per_oct; ++n) {
    // from 0 to 4, 4 to 5, 5 to 1, 1 to 0
    int ang1 = 0 * nang_per_oct + n;
    int ang2 = 3 * nang_per_oct + n;

    int ang3 = 1 * nang_per_oct + n;
    int ang4 = 2 * nang_per_oct + n;

    Real temp = ir[ang1];
    ir[ang1] = ir[ang2];
    ir[ang2] = temp;


    temp = ir[ang3];
    ir[ang3] = ir[ang4];
    ir[ang4] = temp;
  }
  return;
}

// Claude: Rotation by pi about the z axis: (mu_x, mu_y) -> (-mu_x, -mu_y), i.e. octants
// 0 <-> 3, 1 <-> 2 (z > 0) and 4 <-> 7, 5 <-> 6 (z < 0). Both sign flips are mirror
// operations, so angle i of one octant maps to angle i of the other; and a rotation by
// pi is its own inverse, so using the same swap at the inner and outer faces is
// correct. 'direction' is ignored.
void CopyIntensity3(Real *ir, int nang_per_oct, int direction) {
  // here ir is only intensity for each cell and each frequency band
  for (int n=0; n<nang_per_oct; ++n) {
    // switch 1 and 2, switch 3 and 0
    int ang1 = 0 * nang_per_oct + n;
    int ang2 = 3 * nang_per_oct + n;
    int ang3 = 1 * nang_per_oct + n;
    int ang4 = 2 * nang_per_oct + n;
    Real temp = ir[ang1];
    ir[ang1] = ir[ang2];
    ir[ang2] = temp;

    temp = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = temp;
  }
  for (int n=0; n<nang_per_oct; ++n) {
    // switch 5 and 6, switch 4 and 7
    int ang1 = 5 * nang_per_oct + n;
    int ang2 = 6 * nang_per_oct + n;
    int ang3 = 4 * nang_per_oct + n;
    int ang4 = 7 * nang_per_oct + n;
    Real temp = ir[ang1];
    ir[ang1] = ir[ang2];
    ir[ang2] = temp;
    temp = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = temp;
  }
  return;
}

// The angular octant ( in x-y plane) is
//   1  |  0       5  |  4
//   -------      ---------
//   3  |  2       7  |  6

// Claude: Quarter turn about the z axis, as a 4-cycle of octants in each z half: new[0]
// = old[2], new[2] = old[3], new[3] = old[1], new[1] = old[0] (and 4..7 likewise). With
// the octant numbering above, that is the relabelling for a ray turned by +pi/2.
//
// Claude: SUSPECT (angle_flag = 0): 'direction' is ignored, and RotateHPi_InnerX3 and
// RotateHPi_OuterX3 both call this with the same sense. The outer ghost cells hold data
// copied from phi - pi/2 and the inner ones from phi + pi/2, so the two faces need
// opposite rotations; one of them gets the wrong one. Also, a quarter turn swaps |mu_x|
// and |mu_y|, which in general maps angle i of an octant to a DIFFERENT index of the
// target octant, but this keeps i. That is only right if each octant's angles happen to
// be ordered symmetrically under mu_x <-> mu_y.
void CopyIntensity4(Real *ir, int nang_per_oct, int direction) {
  // here ir is only intensity for each cell and each frequency band
  for (int n=0; n<nang_per_oct; ++n) {
    // from 0 to 1, 1 to 3, 3 to 2, 2 to 0
    int ang1 = 0 * nang_per_oct + n;
    int ang2 = 1 * nang_per_oct + n;
    int ang3 = 3 * nang_per_oct + n;
    int ang4 = 2 * nang_per_oct + n;
    Real temp = ir[ang1];
    ir[ang1] = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = ir[ang2];
    ir[ang2] = temp;
  }
  for (int n=0; n<nang_per_oct; ++n) {
    // from 4 to 5, 5 to 7, 7 to 6, 6 to 1
    int ang1 = 4 * nang_per_oct + n;
    int ang2 = 5 * nang_per_oct + n;
    int ang3 = 7 * nang_per_oct + n;
    int ang4 = 6 * nang_per_oct + n;

    Real temp = ir[ang1];
    ir[ang1] = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = ir[ang2];
    ir[ang2] = temp;
  }
  return;
}

// Claude: Used for rotate_theta at the inner x2 face. It swaps every octant with the
// octant of opposite mu_z sign (0 <-> 4, 1 <-> 5, 2 <-> 6, 3 <-> 7), i.e. mu_z -> -mu_z
// with the same index i inside the octant. 'direction' is ignored.
void CopyIntensity5(Real *ir, int nang_per_oct, int direction) {
  // here ir is only intensity for each cell and each frequency band
  for (int n=0; n<nang_per_oct; ++n) {
    // switch 0 and 4, switch 5 and 1
    int ang1 = 0 * nang_per_oct + n;
    int ang2 = 4 * nang_per_oct + n;
    int ang3 = 5 * nang_per_oct + n;
    int ang4 = 1 * nang_per_oct + n;
    Real temp = ir[ang1];
    ir[ang1] = ir[ang2];
    ir[ang2] = temp;

    temp = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = temp;
  }
  for (int n=0; n<nang_per_oct; ++n) {
    // Claude: This loop swaps 2 <-> 6 and 3 <-> 7; the comment on the next line
    // describes a cycle 4 -> 5 -> 7 -> 6 that the code does not do.
    // from 4 to 5, 5 to 7, 7 to 6, 6 to 4
    int ang1 = 2 * nang_per_oct + n;
    int ang2 = 6 * nang_per_oct + n;
    int ang3 = 7 * nang_per_oct + n;
    int ang4 = 3 * nang_per_oct + n;
    Real temp = ir[ang1];
    ir[ang1] = ir[ang2];
    ir[ang2] = temp;

    temp = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = temp;
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x2, inner x2 boundary by Pi/2

// Claude: Called only when rotate_theta == 1. Permutes the angles of every ghost cell
// below jl in place, after the periodic copy has filled them.
// anti-clockwise rotation
// This function is used after periodic copy

void RadBoundaryVariable::RotateHPi_InnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones
  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(k,jl-j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity5(ir, nang_per_oct, 1);
        }
      }
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  ROTATE boundary radiation conditions, outer x2 boundary

// Claude: Not called anywhere: ApplyPhysicalBoundaries and ApplyRadPhysicalBoundaries
// only ever call RotateHPi_InnerX2 for rotate_theta.
void RadBoundaryVariable::RotateHPi_OuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct =pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(k,ju+j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity5(ir, nang_per_oct, -1);
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_InnerX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, inner x3 boundary by Pi/2

// Claude: rotate_phi == 1: quarter-turn relabelling of the ghost cells below kl. See
// the SUSPECT note on CopyIntensity4 about the rotation sense.
// anti-clockwise rotation
// This function is used after periodic copy
void RadBoundaryVariable::RotateHPi_InnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {
  // copy radiation variables into ghost zones
  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq =pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(kl-k,j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity4(ir, nang_per_oct, 1);
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_OuterX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, outer x3 boundary by Pi/2

// Claude: rotate_phi == 1: same CopyIntensity4 call as the inner face, so both faces
// rotate the same way. See the SUSPECT note on CopyIntensity4.
// anti-clockwise rotation
// This function is used after periodic copy
void RadBoundaryVariable::RotateHPi_OuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {
  // copy radiation variables into ghost zones
  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(ku+k,j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity4(ir, nang_per_oct, 1);
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotatePi_InnerX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, inner x3 boundary by Pi

// Claude: rotate_phi == 2: rotation by pi of the ghost cells below kl.
// anti-clockwise rotation
// This function is used after periodic copy

void RadBoundaryVariable::RotatePi_InnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {
  // copy radiation variables into ghost zones
  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(kl-k,j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity3(ir, nang_per_oct, 1);
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotatePi_OuterX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, outer x3 boundary by Pi

// Claude: rotate_phi == 2: rotation by pi of the ghost cells above ku (same swap is
// correct).
// anti-clockwise rotation
// This function is used after periodic copy
void RadBoundaryVariable::RotatePi_OuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {
  // copy radiation variables into ghost zones
  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(ku+k,j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity3(ir, nang_per_oct, 1);
        }
      }
    }
  }
  return;
}
