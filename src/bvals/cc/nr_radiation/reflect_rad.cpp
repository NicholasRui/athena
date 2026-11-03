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


// Claude: Overview. A reflecting boundary mirrors the solution across the boundary
// face. Ghost cell il-i takes the data of interior cell il+i-1 (its mirror image in
// space), and each ray in the ghost cell takes the intensity of the interior ray whose
// direction is the mirror image across the face: the component of mu normal to the face
// flips sign, the other two are unchanged.
//
// Claude: These routines find the mirrored ray by swapping whole OCTANT BLOCKS. That
// only works for the angle_flag = 0 / 10 layout, where angles are stored octant by
// octant, n = octant*nang_per_oct + i, with octants numbered as in the diagram below
// (bit 0 = sign of mu_x, bit 1 = sign of mu_y, bit 2 = sign of mu_z), and where angle i
// of one octant is the mirror image of angle i of the neighbouring octant.
//
// Claude: BUG (angle_flag = 1): AngularGrid(angle_flag, nzeta, npsi) stores angles as n
// = izeta*(2*npsi) + ipsi, with the two polar angles (if any) at nang-2 and nang-1.
// There are no octant blocks, so nang/noct means nothing in that layout and these swaps
// pair up unrelated rays. With polar_angle the integer division nang/noct also drops
// the 2 polar angles, so they are never reflected at all. The correct mirrored index
// for angle_flag = 1 in spherical polar (mu components are local r, theta, phi) is: x1
// (r) face: mu_r flips, so izeta -> 2*nzeta-1-izeta with ipsi unchanged, and the polar
// angles swap (nang-2 <-> nang-1). x2 (theta) face: mu_theta = sin(zeta)*cos(psi)
// flips, so psi -> pi - psi, i.e. ipsi -> (npsi-1-ipsi) mod 2*npsi; except npsi == 1,
// where AngularGrid puts the transverse component on mu_theta itself, so ipsi ->
// 1-ipsi. Polar angles unchanged. x3 (phi) face: mu_phi = sin(zeta)*sin(psi) flips, so
// ipsi -> 2*npsi-1-ipsi. A layout-independent fix is to precompute, once per face, the
// index of the angle whose mu equals the reflected mu, and use that table here.
// The angular grid can change
// The angular octant is
//   1  |  0       5  |  4
//   -------      ---------
//   3  |  2       7  |  6
// in radiation class, nang_per_oct is angles per octant, noct is the number of octant
// radiation relfection means set specific intensitiy corresponding
// the value in the opposite octant

// Claude: For one cell and one frequency group, swap octant blocks li and lo between
// the interior cell (iri) and the ghost cell (iro): ghost rays in octant li get the
// interior rays of octant lo, and vice versa. iri and iro point at the first of this
// group's nang angles, so any index >= nang lands in the NEXT frequency group of the
// same cell, or, for the last group, in the next cell along x1.
// Temporary function to copy intensity
void CopyIntensity(Real *iri, Real *iro, int li, int lo, int nang_per_oct) {
  // here ir is only intensity for each cell and each frequency band
  for (int n=0; n<nang_per_oct; ++n) {
    int angi = li * nang_per_oct + n;
    int ango = lo * nang_per_oct + n;
    iro[angi] = iri[ango];
    iro[ango] = iri[angi];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x1 boundary

void RadBoundaryVariable::ReflectInnerX1(Real time, Real dt, int il, int jl, int ju,
                                         int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  // Claude: nang_per_oct = nang/noct assumes the octant layout described at the top of
  // the file.
  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // Claude: i = 1 is the ghost cell touching the face; it mirrors the first active
      // cell il.
      for (int i=1; i<=ngh; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *iri = &var(k,j,(il+i-1),ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(k,j, il-i, ifr*pmy_block_->pnrrad->nang);
          // Claude: Flip mu_x: swap octants that differ only in the sign of mu_x. 0 <->
          // 1 always; 2 <-> 3 in 2D and 3D; 4 <-> 5 and 6 <-> 7 in 3D only.
          CopyIntensity(iri, iro, 0, 1, nang_per_oct);
          if (noct > 2) {
            CopyIntensity(iri, iro, 2, 3, nang_per_oct);
          }
          // Claude: BUG: should be noct > 4. In 2D noct == 4, so octants 4-7 do not
          // exist. Their indices (>= 4*nang_per_oct = nang) run past this frequency
          // group: into the next group, or for the last group into the next cell. For
          // the ghost cell touching the face (i = 1) that next cell is the ACTIVE cell
          // il, so active intensities get overwritten. Confirmed in a 2D Cartesian
          // angle_flag = 0 run with reflecting x1 boundaries and zero time steps: the
          // first active column changed, and the run aborted at exit with heap
          // corruption (munmap_chunk(): invalid pointer).
          if (noct > 3) {
            CopyIntensity(iri, iro, 4, 5, nang_per_oct);
            CopyIntensity(iri, iro, 6, 7, nang_per_oct);
          }
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX1(
//          Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x1 boundary

void RadBoundaryVariable::ReflectOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          // Claude: Ghost cell iu+i mirrors interior cell iu-i+1.
          Real *iri = &var(k,j,(iu-i+1),ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(k,j, iu+i, ifr*pmy_block_->pnrrad->nang);
          // Claude: Flip mu_x, as for the inner x1 face.
          CopyIntensity(iri, iro, 0, 1, nang_per_oct);
          if (noct > 2) {
            CopyIntensity(iri, iro, 2, 3, nang_per_oct);
          }
          // Claude: BUG: should be noct > 4 (same as ReflectInnerX1). Here the overrun
          // goes past the ghost cell iu+i into iu+i+1; for i = ngh that is the first
          // cell of the next j row, or past the end of the array on the last row.
          if (noct > 3) {
            CopyIntensity(iri, iro, 4, 5, nang_per_oct);
            CopyIntensity(iri, iro, 6, 7, nang_per_oct);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x2 boundary

void RadBoundaryVariable::ReflectInnerX2(Real time, Real dt, int il, int iu, int jl,
                                         int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx

  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          // Claude: Ghost cell jl-j mirrors interior cell jl+j-1.
          Real *iri = &var(k,jl+j-1,i,ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(k,jl-j,i, ifr*pmy_block_->pnrrad->nang);
          // Claude: Flip mu_y: swap octants that differ only in the sign of mu_y. 0 <->
          // 2 and 1 <-> 3 always; 4 <-> 6 and 5 <-> 7 in 3D only.
          CopyIntensity(iri, iro, 0, 2, nang_per_oct);
          CopyIntensity(iri, iro, 1, 3, nang_per_oct);

          // Claude: BUG: should be noct > 4. In 2D (noct == 4) octants 4-7 overrun into
          // the next frequency group or the neighbouring cell, as in ReflectInnerX1.
          if (noct > 3) {
            CopyIntensity(iri, iro, 4, 6, nang_per_oct);
            CopyIntensity(iri, iro, 5, 7, nang_per_oct);
          }
        }
      }
    }
  }
  return;
}



//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x2 boundary

void RadBoundaryVariable::ReflectOuterX2(Real time, Real dt, int il, int iu, int ju,
                                         int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx

  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          // Claude: Ghost cell ju+j mirrors interior cell ju-j+1.
          Real *iri = &var(k,ju-j+1,i,ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(k,ju+j,i, ifr*pmy_block_->pnrrad->nang);
          // Claude: Flip mu_y, as for the inner x2 face.
          CopyIntensity(iri, iro, 0, 2, nang_per_oct);
          CopyIntensity(iri, iro, 1, 3, nang_per_oct);

          // Claude: BUG: should be noct > 4, as in ReflectInnerX2.
          if (noct > 3) {
            CopyIntensity(iri, iro, 4, 6, nang_per_oct);
            CopyIntensity(iri, iro, 5, 7, nang_per_oct);
          }
        }
      }
    }
  }
  return;
}




//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX3(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x3 boundary

void RadBoundaryVariable::ReflectInnerX3(Real time, Real dt, int il, int iu, int jl,
                                         int ju, int kl, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx

  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          // Claude: Ghost cell kl-k mirrors interior cell kl+k-1.
          Real *iri = &var(kl+k-1,j,i,ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(kl-k,j,i, ifr*pmy_block_->pnrrad->nang);
          // Claude: Flip mu_z: swap octants that differ only in the sign of mu_z (0 <->
          // 4, 1 <-> 5, 2 <-> 6, 3 <-> 7). x3 faces only exist in 3D, where noct == 8,
          // so these are in range for angle_flag = 0.
          CopyIntensity(iri, iro, 0, 4, nang_per_oct);
          CopyIntensity(iri, iro, 1, 5, nang_per_oct);
          CopyIntensity(iri, iro, 2, 6, nang_per_oct);
          CopyIntensity(iri, iro, 3, 7, nang_per_oct);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX3(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x3 boundary

void RadBoundaryVariable::ReflectOuterX3(Real time, Real dt, int il, int iu, int jl,
                                         int ju, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  const int& noct = pmy_block_->pnrrad->noct;
  int nang_per_oct = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          // Claude: Ghost cell ku+k mirrors interior cell ku-k+1.
          Real *iri = &var(ku-k+1,j,i,ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(ku+k,j,i, ifr*pmy_block_->pnrrad->nang);
          // Claude: Flip mu_z, as for the inner x3 face.
          CopyIntensity(iri, iro, 0, 4, nang_per_oct);
          CopyIntensity(iri, iro, 1, 5, nang_per_oct);
          CopyIntensity(iri, iro, 2, 6, nang_per_oct);
          CopyIntensity(iri, iro, 3, 7, nang_per_oct);
        }
      }
    }
  }
  return;
}
