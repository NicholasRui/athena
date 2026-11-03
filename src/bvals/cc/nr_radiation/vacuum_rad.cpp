//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file vacuum_rad.cpp
//  \brief implementation of vacuum BCs in each dimension

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../nr_radiation/radiation.hpp"
#include "./bvals_rad.hpp"


// Claude: Overview. A vacuum boundary lets radiation leave the domain but nothing come
// in. In the ghost cells, rays pointing out of the domain copy the boundary cell's
// intensity (free streaming outward), and rays pointing into the domain are set to 0
// (no radiation arrives from outside). The test uses the sign of mu, not octant blocks,
// so it works for any angle layout, including angle_flag = 1. With angle_flag = 1,
// mu(0..2) are the local (r, theta, phi) components in spherical polar, or (x, y, z) in
// Cartesian, so testing mu(0), mu(1), mu(2) at the x1, x2, x3 faces is correct.
//
// Claude: Note: rays exactly parallel to a face (mu component == 0) count as incoming
// and are zeroed. With angle_flag = 1 that includes the two polar angles at x2 and x3
// faces, and every ray at x3 faces when npsi == 1. This is harmless for the fluxes
// through the face, since those rays have zero normal velocity.
//
// Claude: The comments 'angles per octant' on nang below are wrong: nang is the total
// number of angles.
// The angular grid can change
// The angular octant is
//   1  |  0       5  |  4
//   -------      ---------
//   3  |  2       7  |  6
// in radiation class, nang_per_oct is angles per octant, noct is the number of octant
// radiation relfection means set specific intensitiy corresponding
// the value in the opposite octant

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::VacuumInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief VACUUM boundary conditions, inner x1 boundary

void RadBoundaryVariable::VacuumInnerX1(Real time, Real dt, int il, int jl, int ju,
                                        int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            // Claude: At the inner x1 face the outward normal is -x1, so mu_x < 0 means
            // the ray is leaving the domain: copy it. mu_x >= 0 would enter the domain:
            // set it to 0. The sign is taken at the boundary cell il and reused for all
            // ngh ghost cells (mu is the same there for angle_flag = 1; for angle_flag
            // = 0 in curvilinear coordinates it varies a little from cell to cell).
            const Real& miux=pmb->pnrrad->mu(0,k,j,il,n);
            if (miux < 0.0) {
              (*var_cc)(k,j,il-i,ang) = (*var_cc)(k,j,il,ang);
            } else {
              (*var_cc)(k,j,il-i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::VacuumOuterX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief VACUUM boundary conditions, outer x1 boundary

void RadBoundaryVariable::VacuumOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {

  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            // Claude: Outer x1 face: the outward normal is +x1, so mu_x > 0 rays leave
            // and are copied.
            const Real& miux=pmb->pnrrad->mu(0,k,j,iu,n);
            if (miux > 0.0) {
              (*var_cc)(k,j,iu+i,ang) = (*var_cc)(k,j,iu,ang);
            } else {
              (*var_cc)(k,j,iu+i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}





//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::VacuumInnerX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
// Claude: The \brief above says REFLECTING; this is the VACUUM inner x2 boundary.
//  \brief REFLECTING boundary conditions, inner x2 boundary

void RadBoundaryVariable::VacuumInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            // Claude: Inner x2 face: mu_y (mu_theta in spherical polar) < 0 leaves the
            // domain.
            const Real& miuy=pmb->pnrrad->mu(1,k,jl,i,n);
            if (miuy < 0.0) {
              (*var_cc)(k,jl-j,i,ang) = (*var_cc)(k,jl,i,ang);
            } else {
              (*var_cc)(k,jl-j,i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::VacuumOuterX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief VACUUM boundary conditions, outer x2 boundary

void RadBoundaryVariable::VacuumOuterX2(
    Real time, Real dt, int il, int iu, int ju,  int kl, int ku, int ngh) {

  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            // Claude: Outer x2 face: mu_y (mu_theta) > 0 leaves the domain.
            const Real& miuy=pmb->pnrrad->mu(1,k,ju,i,n);
            if (miuy > 0.0) {
              (*var_cc)(k,ju+j,i,ang) = (*var_cc)(k,ju,i,ang);
            } else {
              (*var_cc)(k,ju+j,i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}






//----------------------------------------------------------------------------------------
// Claude: The \fn above says VacuumInnerX2; this is VacuumInnerX3.
//! \fn void RadBoundaryVariable::VacuumInnerX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief VACUUM boundary conditions, inner x2 boundary

void RadBoundaryVariable::VacuumInnerX3(
    Real time, Real dt, int il, int iu, int jl,  int ju, int kl, int ngh) {

  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            // Claude: Inner x3 face: mu_z (mu_phi in spherical polar) < 0 leaves the
            // domain.
            const Real& miuz=pmb->pnrrad->mu(2,kl,j,i,n);
            if (miuz < 0.0) {
              (*var_cc)(kl-k,j,i,ang) = (*var_cc)(kl,j,i,ang);
            } else {
              (*var_cc)(kl-k,j,i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::VacuumOuterX3(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief VACUUM boundary conditions, outer X3 boundary

void RadBoundaryVariable::VacuumOuterX3(
    Real time, Real dt, int il, int iu, int jl,  int ju, int ku, int ngh) {

  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            // Claude: Outer x3 face: mu_z (mu_phi) > 0 leaves the domain.
            const Real& miuz=pmb->pnrrad->mu(2,ku,j,i,n);
            if (miuz > 0.0) {
              (*var_cc)(ku+k,j,i,ang) = (*var_cc)(ku,j,i,ang);
            } else {
              (*var_cc)(ku+k,j,i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}
