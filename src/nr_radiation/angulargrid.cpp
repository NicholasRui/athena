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
//! \file angulargrid.cpp
//  \brief implementation of angular grid in class Radiation
//======================================================================================

// C headers

// C++ headers
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "./radiation.hpp"

//--------------------------------------------------------------------------------------
// \!fn void Radiation::AngularGrid(int angle_flag, int nmu)
// \brief function to create the angular grid

void NRRadiation::AngularGrid(int angle_flag, int nmu) {
  std::stringstream msg;
  MeshBlock *pmb=pmy_block;
  // allocate some temporaray arrays
  AthenaArray<Real> mu2tmp, mutmp, wtmp2, wtmp;
  AthenaArray<Real> pmat, pinv, wpf;
  AthenaArray<int> plab, pl;

  if (polar_angle) {
    msg << "### FATAL ERROR in function [InitialAngle]" << std::endl
        << "Polar angle cannot be used with: "<< angle_flag << "\n ";
        throw std::runtime_error(msg.str().c_str());
  }

  int nang_per_oct = nang/noct;

  mu2tmp.NewAthenaArray(nmu);
  mutmp.NewAthenaArray(nang_per_oct,3);
  wtmp2.NewAthenaArray(nmu-1);
  wtmp.NewAthenaArray(nmu);

  pmat.NewAthenaArray(nmu,nmu);
  pinv.NewAthenaArray(nmu-1,nmu-1);
  plab.NewAthenaArray(nang_per_oct);
  pl.NewAthenaArray(nmu,3);
  wpf.NewAthenaArray(nmu-1);


  // initialize coordinate direction
  int axisx=0, axisy=1, axisz=2;
  pmy_block->pcoord->AxisDirection(&axisx, &axisy, &axisz);

  // check the dimension of the problem
  int n1z = pmb->ncells1, n2z = pmb->ncells2, n3z = pmb->ncells3;

  int ndim = 1;
  if (n2z > 1) ndim = 2;
  if (n3z > 1) ndim = 3;

  if (angle_flag == 0) {
    if (ndim > 1) {
      Real deltamu = 2.0 / (2 * nmu - 1);
      mu2tmp(0) = 1.0 / (3.0 * (2 * nmu - 1));
      for (int i=1; i<nmu; i++)
        mu2tmp(i) = mu2tmp(i-1) + deltamu;

      Real w2 = 4.0 * mu2tmp(0);
      Real wsum2 = std::sqrt(w2);
      wtmp2(0) = wsum2;

      for (int i=1; i<nmu-2; i++) {
        w2 += deltamu;
        wtmp2(i) = std::sqrt(w2);
        wsum2 += wtmp2(i);
      }

      if (nmu > 2)
        wtmp2(nmu-2) = 2.0*(nmu-1)/3.0 - wsum2;

      wtmp(0) = wtmp2(0);
      Real wsum = wtmp(0);
      for (int i=1; i<nmu-1; ++i) {
        wtmp(i) = wtmp2(i) - wtmp2(i-1);
        wsum += wtmp(i);
      }
      wtmp(nmu-1) = 1.0 - wsum;

      int np = 0;
      int iang = 0;

      // initialize to be zero
      for (int i=0; i<nmu; ++i)
        for (int j=0; j<nmu; ++j)
          pmat(i,j) = 0.0;

      for (int i=0; i<nmu; i++) {
        for (int j=0; j<nmu; j++) {
          for (int k=0; k<nmu; k++) {
            if (i + j + k == nmu - 1) {
              // assign cosines to temporary array grid
              mutmp(iang,0) = std::sqrt(mu2tmp(j));
              mutmp(iang,1) = std::sqrt(mu2tmp(k));
              mutmp(iang,2) = std::sqrt(mu2tmp(i));

              int ip=Permutation(i,j,k,np,pl);
              if (ip == -1) {
                pl(np,0) = i;
                pl(np,1) = j;
                pl(np,2) = k;
                pmat(i,np) += 1.0;
                plab(iang) = np;
                np++;
              } else {
                pmat(i,ip) += 1.0;
                plab(iang) = ip;
              }

              iang++;
            }
          }
        }
      }
      if (nmu > 1) {
        //  Invert matrix of Permutations families */
        InverseMatrix(nmu-1, pmat,pinv);
        // Solve for and assign weights for each Permutation family
        MatrixMult(nmu-1,nmu-1,pinv,wtmp,wpf);

        for (int l=0; l<noct; ++l)
          for (int i=0; i<nang_per_oct; ++i)
              wmu(l*nang_per_oct+i) = wpf(plab(i));

      } else {
        for (int l=0; l<noct; ++l)
          wmu(l) = 1.0;
      }
    }

    if (ndim == 1) {
      AthenaArray<Real> mutmp1d, wtmp1d;
      mutmp1d.NewAthenaArray(2*nmu);
      wtmp1d.NewAthenaArray(2*nmu);

      Gauleg(2*nmu, -1.0, 1.0, mutmp1d, wtmp1d);
      for (int i=0; i<2*nmu; ++i)
        wmu(i) = 0.5 * wtmp1d(i);

      for (int n1=0; n1<n1z; ++n1)
        for (int i=0; i<2*nmu; ++i)
          mu(0,0,0,n1,i)=mutmp1d(i);
      
    } else if (ndim == 2) {
      // for spherical coordinate system, it should be r-theta-phi
      for (int n2=0; n2<n2z; ++n2) {
        for (int n1=0; n1<n1z; ++n1) {
          for (int j=0; j<2; ++j) {
            for (int k=0; k<2; ++k) {
              int l=2*j+k;
              for (int i=0; i<nang_per_oct; ++i) {
                int mi = l*nang_per_oct + i;
                if (k == 0) {
                  mu(axisx,0,n2,n1,mi) =  mutmp(i,0);
                } else {
                  mu(axisx,0,n2,n1,mi) = -mutmp(i,0);
                }
                if (j == 0) {
                  mu(axisy,0,n2,n1,mi) =  mutmp(i,1);
                } else {
                  mu(axisy,0,n2,n1,mi) = -mutmp(i,1);
                }
              }
            }
          }
        }
      }
      // for the angular weight
      for (int i=0; i<noct * nang_per_oct; ++i)
        wmu(i) *= 0.25;

    } else if (ndim == 3) {
      for (int n3=0; n3<n3z; ++n3) {
        for (int n2=0; n2<n2z; ++n2) {
          for (int n1=0; n1<n1z; ++n1) {
            for (int j=0; j<2; ++j) {
              for (int k=0; k<2; ++k) {
                for (int l=0; l<2; ++l) {
                  int m=4*j+2*k+l;

                  for (int i=0; i<nang_per_oct; ++i) {
                    int mi = m*nang_per_oct + i;
                    if (l == 0) {
                      mu(axisx,n3,n2,n1,mi) =  mutmp(i,0);
                    } else {
                      mu(axisx,n3,n2,n1,mi) = -mutmp(i,0);
                    }
                    if (k == 0) {
                      mu(axisy,n3,n2,n1,mi) =  mutmp(i,1);
                    } else {
                      mu(axisy,n3,n2,n1,mi) = -mutmp(i,1);
                    }
                    if (j == 0) {
                      mu(axisz,n3,n2,n1,mi) =  mutmp(i,2);
                    } else {
                      mu(axisz,n3,n2,n1,mi) = -mutmp(i,2);
                    }
                  } // end i
                } // end l
              } // end k
            } // end j
          } // end n1
        } // end n2
      } // end n3

      for (int i=0; i<noct * nang_per_oct; ++i)
        wmu(i) *= 0.125;
      
    } // end nDIM = 3
  } else if (angle_flag == 10) {
    if ((ndim == 1) || (ndim == 3)) {
      msg <<"[Warning]: ang_quad = 10 should be used" <<
      "only for 2D problems. \n" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    Real delmu=1.0/std::sqrt((Real)3.0);
    Real sintheta = std::sqrt((Real)2.0)/std::sqrt((Real)3.0);
    Real phi = 0.5 * PI / (Real) (2*nmu);

    for (int i=0; i<nmu; ++i) {
      mutmp(i,0) = sintheta * cos(phi*(Real)(2*i+1));
      mutmp(i,1) = sintheta * sin(phi*(Real)(2*i+1));
      mutmp(i,2) = delmu;
    }
    for (int n2=0; n2<n2z; ++n2) {
      for (int n1=0; n1<n1z; ++n1) {
        for (int j=0; j<2; ++j) {
          for (int k=0; k<2; ++k) {
            int l=2*j+k;

            for (int i=0; i<nang_per_oct; ++i) {
              int mi = l*nang_per_oct + i;
              if (k == 0) {
                mu(axisx,0,n2,n1,mi) =  mutmp(i,0);
              } else {
                mu(axisx,0,n2,n1,mi) = -mutmp(i,0);
              }

              if (j == 0) {
                mu(axisy,0,n2,n1,mi) =  mutmp(i,1);
              } else {
                mu(axisy,0,n2,n1,mi) = -mutmp(i,1);
              }
            }
          }
        }
      }
    }
    for (int i=0; i<noct*nang_per_oct; ++i)
      wmu(i) = 0.25/((Real)nmu);
    
  } else {
    msg << "### FATAL ERROR in function [InitialAngle]" << std::endl
    << "Type of angular discretization unknow: "<< angle_flag << "\n ";
    throw std::runtime_error(msg.str().c_str());
  }

  // Now change the angle cosines to different coordinate systems
  pmy_block->pcoord->ConvertAngle(pmy_block,nang, mu);
  return;
}

// overwrite function to create angles that are defined with respect
// to local coordinate system

// This is the so-called spherical polar tetrad,
// which is ideal for spherical polar coordinates

// This is the angle_flag == 1 case and follows Jiang 2022, ApJS 263, 4

void NRRadiation::AngularGrid(int angle_flag, int nzeta, int npsi) {
  MeshBlock *pmb=pmy_block;

  // require angle_flag = 1
  if (angle_flag != 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in AngularGrid(nzeta, npsi)" << std::endl
        << "This overload only handles angle_flag == 1, got: " << angle_flag
        << std::endl;
    ATHENA_ERROR(msg);
  }

  // initialize coordinate directions
  // IMPORTANT: if spherical polar, AxisDirection in
  // src/coordinates/spherical_polar.cpp rewrites (*axisx, *axisy, *axisz) = (1, 2, 0)
  int axisx=0, axisy=1, axisz=2;
  pmy_block->pcoord->AxisDirection(&axisx, &axisy, &axisz);

  // Get number of cells (including ghost) and dimensions
  int n1z = pmb->ncells1, n2z = pmb->ncells2, n3z = pmb->ncells3;

  int ndim = 1;
  if (n2z > 1) ndim = 2;
  if (n3z > 1) ndim = 3;

  // Define max of zeta grid
  // if polar_angle, exclude poles since we will populate these later on
  Real coszeta_max;
  if (polar_angle) {
    coszeta_max = ((Real)nang-2.0)/(Real)nang;
  } else {
    coszeta_max = 1.0;
  }

  // sin(zeta) at cell centers, used when constructing mu below
  AthenaArray<Real> sinzeta_v;

  // separate ghost zones and active zones
  if (nzeta > 0) { // general case: gridding over zeta
    // =======================
    // ====== ZETA GRID ======
    // =======================

    // Compute equally-spaced FACE-CENTERED grid in cos(zeta) [angular grid in theta]
    // write some distinguished values and assert symmetry upon reflection
    // across equator to avoid roundoff-error-induced asymmetries
    coszeta_f.NewAthenaArray(2*nzeta+1);
    dcoszeta_f.NewAthenaArray(2*nzeta);

    coszeta_f(0) = coszeta_max;        // north pole
    coszeta_f(nzeta) = 0.0;            // equator
    coszeta_f(2*nzeta) = -coszeta_max; // south pole

    Real dcoszeta = coszeta_max/nzeta;

    for (int i=1; i<nzeta; ++i) {
      coszeta_f(i) = coszeta_max - i*dcoszeta; // upper hemisphere
      coszeta_f(2*nzeta-i) = -coszeta_f(i);    // lower hemisphere
    }

    // save face-centered grid spacings
    for (int i=0; i<2*nzeta; ++i)
      dcoszeta_f(i) = coszeta_f(i) - coszeta_f(i+1);


    // Compute CELL-CENTERED cos(zeta) grid by averaging face-centered values

    // if polar_angle is not enabled, the discrete <cos^2(zeta)> != 1/3 exactly
    // apply a normalization to force it to be 1/3, to match isotropic limit
    Real norm = polar_angle ? 1.0 : 2*nzeta/std::sqrt(4*nzeta*nzeta-1);

    coszeta_v.NewAthenaArray(2*nzeta);

    for (int i=0; i<nzeta; ++i) {
      coszeta_v(i) = 0.5*(coszeta_f(i)+coszeta_f(i+1))*norm; // upper hemisphere
      coszeta_v(2*nzeta-i-1) = -coszeta_v(i);                // lower hemisphere
    }

    // also save sin(zeta) for convenience later
    sinzeta_v.NewAthenaArray(2*nzeta);
    for (int n=0; n<2*nzeta; ++n)
      sinzeta_v(n) = std::sqrt(1.0 - coszeta_v(n)*coszeta_v(n));

    // Write CELL-CENTERED zeta values INCLUDING GHOST CELLS
    // ghost zone values have same delta(zeta) as active cells closest to boundary
    zeta_v_full.NewAthenaArray(2*nzeta+2*NGHOST);
    dzeta_v_full.NewAthenaArray(2*nzeta+2*NGHOST);

    for (int i=0; i<2*nzeta; ++i) // active cells
      zeta_v_full(i+NGHOST) = acos(coszeta_v(i));
    for (int i=1; i<= NGHOST; ++i) { // ghost cells
      // below lower boundary
      Real delta_lo = zeta_v_full(i+NGHOST) - zeta_v_full(NGHOST+i-1);
      zeta_v_full(NGHOST-i) = zeta_v_full(NGHOST-i+1) - delta_lo;

      // above upper boundary
      Real delta_hi = zeta_v_full(2*nzeta+NGHOST-i) - zeta_v_full(2*nzeta+NGHOST-i-1);
      zeta_v_full(2*nzeta+NGHOST+i-1) = zeta_v_full(2*nzeta+NGHOST+i-2) + delta_hi;
    }

    // Compute delta(zeta) between neighboring cell-centered cells including ghost cells
    for (int i=0; i<2*nzeta+2*NGHOST-1; ++i)
      dzeta_v_full(i) = zeta_v_full(i+1) - zeta_v_full(i);


    // Write FACE-CENTERED zeta values INCLUDING GHOST CELLS
    // ghost zone values have same delta(zeta) as active cells closest to boundary
    zeta_f_full.NewAthenaArray(2*nzeta+1+2*NGHOST);
    dzeta_f_full.NewAthenaArray(2*nzeta+1+2*NGHOST);

    for (int i=0; i<2*nzeta+1; ++i) // active cells
      zeta_f_full(i+NGHOST) = acos(coszeta_f(i));
    for (int i=1; i<=NGHOST; ++i) { // ghost zones
      // below lower boundary
      Real delta_lo = zeta_f_full(i+NGHOST) - zeta_f_full(NGHOST+i-1);
      zeta_f_full(NGHOST-i) = zeta_f_full(NGHOST-i+1) - delta_lo;

      // above upper boundary
      Real delta_hi = zeta_f_full(2*nzeta+NGHOST-i+1) - zeta_f_full(2*nzeta+NGHOST-i);
      zeta_f_full(2*nzeta+NGHOST+i) = zeta_f_full(2*nzeta+NGHOST+i-1) + delta_hi;
    }

    // Compute delta(zeta) between neighboring face-centered cells including ghost cells
    for (int i=0; i<2*nzeta+2*NGHOST; ++i)
      dzeta_f_full(i) = zeta_f_full(i+1)-zeta_f_full(i);
  } else {
    // (placeholder arrays) Special case: No gridding over zeta (only phi grid)
    coszeta_v.NewAthenaArray(1);    coszeta_v(0) = coszeta_max;
    zeta_v_full.NewAthenaArray(1);  zeta_v_full(0) = 0.0;
    zeta_f_full.NewAthenaArray(1);  zeta_f_full(0) = 0.0;
    dzeta_v_full.NewAthenaArray(1); dzeta_v_full(0) = coszeta_max;
    dzeta_f_full.NewAthenaArray(1); dzeta_f_full(0) = coszeta_max;
    coszeta_f.NewAthenaArray(1);    coszeta_f(0) = coszeta_max;
    dcoszeta_f.NewAthenaArray(1);   dcoszeta_f(0) = coszeta_max;
  }


  // ======================
  // ====== PSI GRID ======
  // ======================

  if (npsi > 0) { // for nontrivial psi grid
    // Compute equally-spaced CELL-CENTERED psi grid (0 to 2 pi)
    // force symmetry exactly to avoid asymmetries due to roundoff error
    psi_v.NewAthenaArray(2*npsi);

    Real dpsi = PI / npsi;

    for (int i=0; i<npsi; ++i) {
      psi_v(i) = (i+0.5)*dpsi; // psi from [0, PI)
      psi_v(2*npsi-1-i) = 2.0*PI - psi_v(i); // [PI, 2 PI)
    }


    // Compute FACE-CENTERED psi grid
    // write some distinguished values and assert symmetry upon reflection
    // across psi = PI to avoid roundoff-error-induced asymmetries
    psi_f.NewAthenaArray(2*npsi+1);
    dpsi_f.NewAthenaArray(2*npsi);

    psi_f(0) = 0.0;
    psi_f(npsi) = PI;
    psi_f(2*npsi) = 2.0*PI;
    for (int i=1; i<npsi; ++i) {
      psi_f(i) = i*dpsi;                   // [0, PI)
      psi_f(2*npsi-i) = 2.0*PI - psi_f(i); // (PI, 2 PI]
    }

    // save face-centered grid spacings
    for (int i=0; i<2*npsi; ++i)
      dpsi_f(i) = psi_f(i+1) - psi_f(i);

    // compute face-centered sin(psi)
    // overwrite some distinguished values to avoid roundoff error
    sinpsi_f.NewAthenaArray(2*npsi+1);

    for (int i=0; i<2*npsi+1; ++i)
      sinpsi_f(i) = sin(psi_f(i));
    sinpsi_f(0)      = 0.0;
    sinpsi_f(npsi)   = 0.0;
    sinpsi_f(2*npsi) = 0.0;


    // Write CELL-CENTERED psi values INCLUDING GHOST CELLS
    // ghost zone values have same delta(psi) as active cells closest to boundary
    psi_v_full.NewAthenaArray(2*npsi+2*NGHOST);
    dpsi_v_full.NewAthenaArray(2*npsi+2*NGHOST);

    for (int i=0; i<2*npsi; ++i) // active cells
      psi_v_full(i+NGHOST) = psi_v(i);
    for (int i=1; i<=NGHOST; ++i) { // ghost cells (periodic -- wraps around)
      psi_v_full(NGHOST-i) = psi_v(2*npsi-i) - 2.0*PI;
      psi_v_full(2*npsi+NGHOST+i-1) = psi_v(i-1) + 2.0*PI;
    }

    // Compute delta(psi) between neighboring cell-centered cells including ghost cells
    for (int i=0; i<2*npsi+2*NGHOST-1; ++i)
      dpsi_v_full(i) = psi_v_full(i+1) - psi_v_full(i);


    // Write FACE-CENTERED psi values INCLUDING GHOST CELLS
    // ghost zone values have same delta(psi) as active cells closest to boundary
    psi_f_full.NewAthenaArray(2*npsi+1+2*NGHOST);
    dpsi_f_full.NewAthenaArray(2*npsi+2*NGHOST+1);

    for (int i=0; i<2*npsi+1; ++i) // active cells
      psi_f_full(i+NGHOST) = psi_f(i);
    for (int i=1; i<=NGHOST; ++i) { // ghost cells (periodic -- wraps around)
      psi_f_full(NGHOST-i) = psi_f(2*npsi-i)-2.0*PI;
      psi_f_full(2*npsi+NGHOST+i) = psi_f(i) + 2.0*PI;
    }

    // Compute delta(psi) between neighboring face-centered cells including ghost cells
    for (int i=0; i<2*npsi+2*NGHOST; ++i)
      dpsi_f_full(i) = psi_f_full(i+1)-psi_f_full(i);
  }

  // =========================================
  // ====== CONSTRUCT FULL ANGULAR GRID ======
  // =========================================

  // mu(x,k,j,i,n) : DIRECTIONAL cosines of angles, cosine wrt. positive coordinate axes
  //  x: axis direction (e.g., axisz)
  //  (k,j,i): spatial indices
  //  n: angle index -- n = izeta*(2*npsi) + ipsi

  bool spherical_polar = (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0);

  // 1D case: zeta grid only
  if (ndim == 1) {
    int x = 0; // only x direction allowed
    int k = 0, j = 0; // all yz indices = 0

    for (int i=0; i<n1z; ++i) { // cells
      for (int n=0; n<2*nzeta; ++n)
        mu(x,k,j,i,n) = coszeta_v(n);

      // polar angles
      if (polar_angle) {
        mu(x,k,j,i,2*nzeta) = 1.0;
        mu(x,k,j,i,2*nzeta+1) = -1.0;
      }
    }


  // 2D special case: one psi point per half, xz plane
  //  angles are: (+sinzeta_v, 0, coszeta) and (-sinzeta, 0, coszeta) for each zeta
  } else if (ndim == 2 && npsi == 1) {
    int k = 0; // z index = 0

    for (int j=0; j<n2z; ++j) { // cells
      for (int i=0; i<n1z; ++i) {
        for (int n=0; n<2*nzeta; ++n) { // angles
          // +sin(zeta) angle
          mu(axisz,k,j,i,2*npsi*n) = coszeta_v(n);
          mu(axisx,k,j,i,2*npsi*n) = sinzeta_v(n);

          // -sin(zeta) angle
          mu(axisz,k,j,i,2*npsi*n+1) = coszeta_v(n);
          mu(axisx,k,j,i,2*npsi*n+1) = -sinzeta_v(n);
        }

        // polar angle
        if (polar_angle) {
          mu(axisz,k,j,i,nang-2) = 1.0;
          mu(axisz,k,j,i,nang-1) = -1.0;
          mu(axisx,k,j,i,nang-2) = 0.0;
          mu(axisx,k,j,i,nang-1) = 0.0;
        }
      }
    }


  // 2D spherical polar case: still need 3D angular grid
  } else if (ndim == 2 && npsi > 1 && spherical_polar) {
    int k = 0; // z index = 0

    for (int j=0; j<n2z; ++j) { // cells
      for (int i=0; i<n1z; ++i) {
        for (int n=0; n<2*nzeta; ++n) { // angles
          for (int m=0; m<2*npsi; ++m) {
            mu(axisx,k,j,i,2*npsi*n+m) = sinzeta_v(n) * cos(psi_v(m));
            mu(axisy,k,j,i,2*npsi*n+m) = sinzeta_v(n) * sin(psi_v(m));
            mu(axisz,k,j,i,2*npsi*n+m) = coszeta_v(n);
          }
        }

        // polar angle
        if (polar_angle) {
          mu(axisz,k,j,i,nang-2) = 1.0;
          mu(axisz,k,j,i,nang-1) = -1.0;
          mu(axisx,k,j,i,nang-2) = 0.0;
          mu(axisx,k,j,i,nang-1) = 0.0;
          mu(axisy,k,j,i,nang-2) = 0.0;
          mu(axisy,k,j,i,nang-1) = 0.0;
        }
      }
    }


  // 2D non-spherical-polar case: angular grid is 2D and restricted to xy plane
  } else if (ndim == 2 && npsi > 1 && !spherical_polar) {
    int k = 0; // z index = 0

    for (int j=0; j<n2z; ++j) { // cells
      for (int i=0; i<n1z; ++i) {
        for (int m=0; m<2*npsi; ++m) { // angles
          mu(axisx,k,j,i,m) = cos(psi_v(m));
          mu(axisy,k,j,i,m) = sin(psi_v(m));
        }

        // polar angle
        if (polar_angle) {
          // +x pole
          mu(axisx,k,j,i,nang-2) = 1.0;
          mu(axisy,k,j,i,nang-2) = 0.0;

          // -x pole
          mu(axisx,k,j,i,nang-1) = -1.0;
          mu(axisy,k,j,i,nang-1) = 0.0;
        }
      }
    }

  // 3D case: full 3D angular grid
  } else if (ndim == 3) {
    for (int k=0; k<n3z; ++k) { // cells
      for (int j=0; j<n2z; ++j) {
        for (int i=0; i<n1z; ++i) {
          for (int n=0; n<2*nzeta; ++n) { // angles
            for (int m=0; m<2*npsi; ++m) {
              mu(axisx,k,j,i,2*npsi*n+m) = sinzeta_v(n) * cos(psi_v(m));
              mu(axisy,k,j,i,2*npsi*n+m) = sinzeta_v(n) * sin(psi_v(m));
              mu(axisz,k,j,i,2*npsi*n+m) = coszeta_v(n);
            }
          }

          // polar angle
          if (polar_angle) {
            mu(axisz,k,j,i,nang-2) = 1.0;
            mu(axisz,k,j,i,nang-1) = -1.0;
            mu(axisx,k,j,i,nang-2) = 0.0;
            mu(axisx,k,j,i,nang-1) = 0.0;
            mu(axisy,k,j,i,nang-2) = 0.0;
            mu(axisy,k,j,i,nang-1) = 0.0;
          }
        }
      }
    }
  }

  // assign angles equal weight (equal subtended solid angle) such that sum weights = 1
  for (int n=0; n<nang; ++n)
    wmu(n) = 1.0/nang;
}
