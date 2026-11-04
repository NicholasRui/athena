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
//! \file radiation.cpp
//  \brief implementation of functions in class Radiation
//======================================================================================

// C headers

// C++ headers
#include <cmath>   // pow
#include <cstdio>  // fopen and fwrite
#include <cstring> // strcmp
#include <iostream>  // cout
#include <sstream>  // msg
#include <stdexcept> // runtime erro
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "implicit/radiation_implicit.hpp"
#include "integrators/rad_integrators.hpp"
#include "radiation.hpp"
// constructor, initializes data structures and parameters

// The default opacity function.
// Do nothing. Keep the opacity as the initial value
inline void DefaultFrequency(NRRadiation *prad) {
  return;
}

// The default opacity function.
// Do nothing. Keep the opacity as the initial value
inline void DefaultEmission(NRRadiation *prad, Real tgas) {
  const int &nfreq = prad->nfreq;
  if (nfreq > 1) {
    for(int ifr=0; ifr<nfreq-1; ++ifr) {
      prad->emission_spec(ifr) =
          prad->IntPlanckFunc(prad->nu_grid(ifr)/tgas, prad->nu_grid(ifr+1)/tgas);
    }
    prad->emission_spec(nfreq-1) = 1.0 - prad->FitIntPlanckFunc(
        prad->nu_grid(nfreq-1)/tgas);
  }
  return;
}


// The default opacity function.
// Do nothing. Keep the opacity as the initial value
inline void DefaultOpacity(MeshBlock *pmb, AthenaArray<Real> &prim) {
  return;
}


NRRadiation::NRRadiation(MeshBlock *pmb, ParameterInput *pin):
    pmy_block(pmb), ir(pmb->ncells3,pmb->ncells2,pmb->ncells1,pmb->nfre_ang),
    ir1(pmb->ncells3,pmb->ncells2,pmb->ncells1,pmb->nfre_ang),
    // constructor overload resolution of non-aggregate class type AthenaArray<Real>
    flux{ {pmb->ncells3, pmb->ncells2, pmb->ncells1+1, pmb->nfre_ang},
          {pmb->ncells3, pmb->ncells2+1, pmb->ncells1, pmb->nfre_ang,
           (pmb->pmy_mesh->f2 ? AthenaArray<Real>::DataStatus::allocated :
            AthenaArray<Real>::DataStatus::empty)},
          {pmb->ncells3+1, pmb->ncells2, pmb->ncells1, pmb->nfre_ang,
           (pmb->pmy_mesh->f3 ? AthenaArray<Real>::DataStatus::allocated :
            AthenaArray<Real>::DataStatus::empty)}
    },
    coarse_ir_(pmb->ncc3, pmb->ncc2, pmb->ncc1,pmb->nfre_ang,
               (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
                AthenaArray<Real>::DataStatus::empty)),
    rad_bvar(pmb, &ir, &coarse_ir_, flux) {
  // universal constants we need
  // https://physics.info/constants/
  // arad = 4 * sigma/c
  Real arad = 7.565733250033928e-15;
  Real c_speed = 2.99792458e10;  // speed of light

  // Read number of cells in each dimension
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;

  int ndim = 1;
  if (nc2 > 1) ndim = 2;
  if (nc3 > 1) ndim = 3;

  bool use_spherical_polar = (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0);

  // read in the parameters
  // angle grid: angle_grid_type (or deprecated angle_flag), and
  // nmu (lab_davis, lab_magic) or nzeta, npsi, polar_angle (rotating)
  // nzeta: total number of polar angles covering 0 to pi/2
  // npsi: total number of azimuthal angles covering 0 to pi
  ReadAngleGridParameters(pin, ndim, angle_flag, nmu, nzeta, npsi, polar_angle);

  rotate_theta=pin->GetOrAddInteger("radiation","rotate_theta",0);
  rotate_phi=pin->GetOrAddInteger("radiation","rotate_phi",0);

  vmax = pin->GetOrAddReal("radiation","vmax",0.9);
  user_unit_ = pin->GetOrAddInteger("radiation","unit",0);
  tunit = pin->GetOrAddReal("radiation","T_unit",1.e7);
  rhounit = pin->GetOrAddReal("radiation","density_unit",1.0);
  lunit = pin->GetOrAddReal("radiation","length_unit",1.0);
  mol_weight = pin->GetOrAddReal("radiation","molecular_weight",0.6);
  kappa_es = pin->GetOrAddReal("radiation","electron_scattering",0.34);

  kappa_es = kappa_es * rhounit * lunit; // dimensionless electron scattering

  if (user_unit_ == 0) {
    prat = pin->GetReal("radiation","prat");
    crat = pin->GetReal("radiation","crat");
  } else if (user_unit_ == 1) {
    // calculate prat and crat based on user provided unit
    Real r_ideal = 8.314462618e7/mol_weight;
    prat = arad * tunit * tunit * tunit/(rhounit * r_ideal);
    Real cs_iso = std::sqrt(r_ideal * tunit);
    crat = c_speed/cs_iso;
  }

  // equivalent temperature for electron
  telectron = 5.94065e9;
  telectron /= tunit;

  reduced_c  = crat * pin->GetOrAddReal("radiation","reduced_factor",1.0);

  Real taucell = pin->GetOrAddReal("radiation","taucell",5);

  //Mesh *pm = pmb->pmy_mesh;
  //ir_output=pin->GetOrAddInteger("radiation","ir_output",0);

  set_source_flag = pin->GetOrAddInteger("radiation","source_flag",1);

  // Assign number of angles (and octants)
  if (angle_flag == ANGLE_FLAG_COROTATING) {
    GetNumAnglesRotating(angle_flag, ndim, polar_angle, nzeta, npsi,
        COORDINATE_SYSTEM, nang);

    // the rotating grid has no octants; noct is only kept for the octant-based
    // reflect/rotate boundaries (reflect_rad.cpp, rotate_rad.cpp), which use nang/noct
    noct = std::pow(2, ndim);
    if (ndim == 2 && use_spherical_polar) {
      noct = 8; // exception for sph coords
    }
  } else {
    GetNumAnglesAndOctantsLab(angle_flag, ndim, nmu, noct, nang);
  }

  // frequency grid
  //frequency grid covers -infty to infty, default nfreq=1, means gray
  // integrated over all frequency

  // the number of frequeny bins
  nfreq  = pin->GetOrAddInteger("radiation","n_frequency",1);
  // when the emission spectrum depends on tgas, we perform

  // minimum frequency
  nu_min = pin->GetOrAddReal("radiation","frequency_min",0);

  // maximum frequency
  nu_max = pin->GetOrAddReal("radiation","frequency_max",1);

  // log frequency space ratio
  fre_ratio= pin->GetOrAddReal("radiation","frequency_ratio",1.0);

  log_fre_ = pin->GetOrAddInteger("radiation","log_frequency_spacing",1);

  restart_from_gray  = pin->GetOrAddInteger("radiation","res_gray",0);
  // when the emission spectrum depends on tgas, we perform

  // setup the frequency grid
  FrequencyGrid();

  UserFrequency = DefaultFrequency;
  UserEmissionSpec = DefaultEmission;

  n_fre_ang = nang * nfreq;
  //co-moving frame frequency grid depends on angels

  // do not add radiation to vars_cc, which needs to be done in different order for
  // restriction/prolongation in AMR
  //  pmb->RegisterMeshBlockData(ir);

  // If user-requested time integrator is type 3S*, allocate additional memory registers
  std::string integrator = pin->GetOrAddString("time", "integrator", "vl2");
  if (integrator == "ssprk5_4" || STS_ENABLED) {
    // future extension may add "int nregister" to Hydro class
    ir2.NewAthenaArray(nc3, nc2, nc1, n_fre_ang);
  }

  ir_old.NewAthenaArray(nc3,nc2,nc1,n_fre_ang);

  if (restart_from_gray) {
    ir_gray.NewAthenaArray(nc3,nc2,nc1,nang);
  }

  // Do not add to cell-centered refinement, as
  // radiation variables need to be done in different order
  //  if (pm->multilevel) {
  //    refinement_idx = pmy_block->pmr->AddToRefinement(&ir, &coarse_ir_);
  //  }

  rad_mom.NewAthenaArray(13,nc3,nc2,nc1);
  rad_mom_cm.NewAthenaArray(4,nc3,nc2,nc1);
  // dump the moments in each frequency groups
  if (nfreq > 1) {
    // moments in different frequency bins
    rad_mom_nu.NewAthenaArray(13*nfreq,nc3,nc2,nc1);
    rad_mom_cm_nu.NewAthenaArray(4*nfreq,nc3,nc2,nc1);
  }

  // the equation is
  // (sigma_s+sigma_a)(J-I)  // Rosseland mean
  // + sigma_p * a_rT^4 - sigma_pe * J // Planck mean

  sigma_s.NewAthenaArray(nc3,nc2,nc1,nfreq);
  sigma_a.NewAthenaArray(nc3,nc2,nc1,nfreq);
  sigma_p.NewAthenaArray(nc3,nc2,nc1,nfreq);
  sigma_pe.NewAthenaArray(nc3,nc2,nc1,nfreq);

  t_floor_.NewAthenaArray(nc3,nc2,nc1);
  t_ceiling_.NewAthenaArray(nc3,nc2,nc1);

  output_sigma.NewAthenaArray(3*nfreq,nc3,nc2,nc1);

  mu.NewAthenaArray(3,nc3,nc2,nc1,nang);
  wmu.NewAthenaArray(nang);

  if (angle_flag == ANGLE_FLAG_COROTATING) {
    AngularGrid(angle_flag, nzeta, npsi);
    if (nc2 > 1) {
      cot_theta.NewAthenaArray(nc2);
      for(int i=0; i<nc2; ++i) {
        Real th_l = pmb->pcoord->x2f(i);
        Real th_r = pmb->pcoord->x2f(i+1);
        cot_theta(i) = (sin(th_r) - sin(th_l)) / (cos(th_l) - cos(th_r));
      }
    }
  } else {
    AngularGrid(angle_flag, nmu);
  }

  // set a default opacity function
  UpdateOpacity = DefaultOpacity;

  pradintegrator = new RadIntegrator(this, pin);

  rad_bvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&rad_bvar);
  // enroll radiation boundary value object
  if (NR_RADIATION_ENABLED) {
    pmb->pbval->bvars_main_int.push_back(&rad_bvar);
  }

  //-----------------------------------------
  // temporary arrays for multi-group moments
  cosx_cm_.NewAthenaArray(nang);
  cosy_cm_.NewAthenaArray(nang);
  cosz_cm_.NewAthenaArray(nang);
  //------------------------------------------

  // set the default t_floor and t_ceiling
  for(int k=0; k<nc3; ++k)
    for(int j=0; j<nc2; ++j)
      for(int i=0; i<nc1; ++i) {
        t_floor_(k,j,i) = TINY_NUMBER;
        t_ceiling_(k,j,i) = HUGE_NUMBER;
      }

  // dump the angular grid and radiation parameters in a file
  if (Globals::my_rank ==0) {
    FILE *pfile;
    std::stringstream msg;
    if ((pfile = fopen("Rad_angles.txt","w")) == NULL) {
      msg << "### FATAL ERROR in Radiation Class" << std::endl
          << "Output file Rad_angles.txt could not be opened";
      throw std::runtime_error(msg.str().c_str());
    }
    // damp the angular grid in one cell

    fprintf(pfile,"prat          %4.2e \n",prat);
    fprintf(pfile,"crat          %4.2e \n",crat);
    fprintf(pfile,"reduced_c     %4.2e \n",reduced_c);
    fprintf(pfile,"Vmax          %4.2e \n",vmax);
    fprintf(pfile,"Tunit         %4.2e \n",tunit);
    fprintf(pfile,"Compt         %d  \n",pradintegrator->compton_flag_);
    fprintf(pfile,"rotate_theta  %d  \n",rotate_theta);
    fprintf(pfile,"rotate_phi    %d  \n",rotate_phi);
    fprintf(pfile,"adv_flag:     %d  \n",pradintegrator->adv_flag_);
    fprintf(pfile,"nzeta:        %d  \n",nzeta);
    fprintf(pfile,"npsi:         %d  \n",npsi);
    fprintf(pfile,"taucell:      %e  \n",taucell);
    fprintf(pfile,"source_flag:  %d  \n",set_source_flag);
    fprintf(pfile,"nfreq      :  %d  \n",nfreq);
    fprintf(pfile,"fre_ratio:    %e  \n",fre_ratio);
    fprintf(pfile,"nu_min:       %e  \n",nu_min);
    fprintf(pfile,"nu_max:       %e  \n",nu_max);
    if (IM_RADIATION_ENABLED) {
      fprintf(pfile,"iteration:    %d  \n",pmb->pmy_mesh->pimrad->ite_scheme);
      fprintf(pfile,"red_or_black: %d  \n",pmb->pmy_mesh->pimrad->rb_or_not);
      fprintf(pfile,"err_limit:    %e  \n",pmb->pmy_mesh->pimrad->error_limit_);
      fprintf(pfile,"n_limit:      %d  \n",pmb->pmy_mesh->pimrad->nlimit_);
      fprintf(pfile,"tau_scheme    %d  \n",pradintegrator->tau_flag_);
    }

    for(int n=0; n<nang; ++n) {
      fprintf(pfile,"%2d   %e   %e   %e    %e\n",n,mu(0,0,0,0,n),mu(1,0,0,0,n),
              mu(2,0,0,0,n), wmu(n));
    }
    if (nfreq > 1) {
      fprintf(pfile,"fre   spec\n");
      for(int ifr=0; ifr<nfreq; ++ifr) {
        fprintf(pfile,"%e   %e\n",nu_grid(ifr),emission_spec(ifr));
      }
    }
    fclose(pfile);
  }
}

NRRadiation::~NRRadiation() {
  delete pradintegrator;
}


void NRRadiation::EnrollOpacityFunction(OpacityFunc MyOpacityFunction) {
  UpdateOpacity = MyOpacityFunction;
}


void NRRadiation::EnrollFrequencyFunction(FrequencyFunc MyFrequencyFunction) {
  UserFrequency = MyFrequencyFunction;
}


void NRRadiation::EnrollEmissionFunction(EmissionFunc MyEmissionSpec) {
  UserEmissionSpec = MyEmissionSpec;
}

// Reads <radiation> parameters for the angular grid
// Called both here and in src/mesh/meshblock.cpp
// Only parameters used by the chosen angle grid are added to pin (GetOrAdd), so that
// later calls (other MeshBlocks, restarts) do not see unused parameters as specified
//
// angle_grid_type sets angle_flag:
// - "lab_davis"  - lab-frame grid of Davis et al. (2012, ApJS, 199, 9)   (angle_flag = 0)
// - "rotating"   - corotating angle grid of Jiang (2022, ApJS, 263, 4)   (angle_flag = 1)
// - "lab_magic"  - lab-frame grid using magic angles for zeta           (angle_flag = 10)
//                  (coszeta = 1/sqrt(3)) for approximate isotropy
// specifying angle_flag manually is still supported, but deprecated and gives warning
//
// lab_davis, lab_magic: nmu required; nzeta, npsi, polar_angle not allowed
// rotating:             nzeta required; npsi required in 2D/3D, not allowed in 1D;
//                       polar_angle optional (0 or 1); nmu not allowed
// parameters not used by the angle grid are returned as 0
void NRRadiation::ReadAngleGridParameters(ParameterInput *pin, int ndim,
    int &angle_flag, int &nmu, int &nzeta, int &npsi, int &polar_angle) {
  std::stringstream msg;
  bool use_spherical_polar = (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0);

  // Assign angle_flag
  if (pin->DoesParameterExist("radiation", "angle_flag")) {
    if (pin->DoesParameterExist("radiation", "angle_grid_type")) {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "Specify only at most one of angle_grid_type and angle_flag (deprecated)"
          << std::endl;
      ATHENA_ERROR(msg);
    }

    angle_flag = pin->GetInteger("radiation", "angle_flag");
    if (angle_flag != ANGLE_FLAG_LAB_DAVIS && angle_flag != ANGLE_FLAG_COROTATING
        && angle_flag != ANGLE_FLAG_LAB_MAGIC) {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "angle_flag = " << angle_flag << " is not one of 0, 1, 10" << std::endl;
      ATHENA_ERROR(msg);
    }

    // warn once per run, not once per MeshBlock
    static bool warned = false;
    if (!warned && Globals::my_rank == 0) {
      std::cout << "### WARNING in radiation class" << std::endl
                << "angle_flag is deprecated; use angle_grid_type = lab_davis (0), "
                << "rotating (1), or lab_magic (10) instead" << std::endl;
    }
    warned = true;

  } else {
    std::string type = pin->GetOrAddString("radiation", "angle_grid_type", "lab_davis");
    if (type == "lab_davis") {
      angle_flag = ANGLE_FLAG_LAB_DAVIS;
    } else if (type == "rotating") {
      angle_flag = ANGLE_FLAG_COROTATING;
    } else if (type == "lab_magic") {
      angle_flag = ANGLE_FLAG_LAB_MAGIC;
    } else {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "angle_grid_type = " << type << " is not one of lab_davis, rotating, "
          << "lab_magic" << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  // Assign nmu, nzeta, npsi, polar_angle
  // GetInteger gives a fatal error if a required parameter is missing
  nmu = 0;
  nzeta = 0;
  npsi = 0;
  polar_angle = 0;

  if (angle_flag == ANGLE_FLAG_LAB_DAVIS || angle_flag == ANGLE_FLAG_LAB_MAGIC) {
    if (pin->DoesParameterExist("radiation", "nzeta")
        || pin->DoesParameterExist("radiation", "npsi")
        || pin->DoesParameterExist("radiation", "polar_angle")) {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "nzeta, npsi, and polar_angle are only used by angle_grid_type = rotating"
          << std::endl;
      ATHENA_ERROR(msg);
    }

    nmu = pin->GetInteger("radiation", "nmu");

    if (nmu < 1) {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "nmu must be >= 1" << std::endl;
      ATHENA_ERROR(msg);
    }

  } else if (angle_flag == ANGLE_FLAG_COROTATING) {
    if (pin->DoesParameterExist("radiation", "nmu")) {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "nmu is not used by angle_grid_type = rotating; use nzeta and npsi"
          << std::endl;
      ATHENA_ERROR(msg);
    }
    // cannot vary angles in psi with only one dimension
    if (ndim == 1 && pin->DoesParameterExist("radiation", "npsi")) {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "npsi is not used by angle_grid_type = rotating in 1D" << std::endl;
      ATHENA_ERROR(msg);
    }

    nzeta = pin->GetInteger("radiation", "nzeta");
    if (ndim > 1) npsi = pin->GetInteger("radiation", "npsi");
    polar_angle = pin->GetOrAddInteger("radiation", "polar_angle", 0);

    if (nzeta < 0 || npsi < 0) {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "nzeta and npsi must be non-negative" << std::endl;
      ATHENA_ERROR(msg);
    }
    if (ndim == 1 && nzeta == 0) {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "nzeta must be >= 1 for angle_grid_type = rotating in 1D" << std::endl;
      ATHENA_ERROR(msg);
    }
    if ((ndim == 3 || (ndim == 2 && use_spherical_polar)) && (nzeta == 0 || npsi == 0)) {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "nzeta and npsi must both be >= 1 for angle_grid_type = rotating "
          << "in 3D or 2D spherical_polar" << std::endl;
      ATHENA_ERROR(msg);
    }
    if (ndim == 2 && nzeta == 0 && npsi == 0) {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "at least one of nzeta and npsi must be >= 1 for angle_grid_type = "
          << "rotating" << std::endl;
      ATHENA_ERROR(msg);
    }
    if (polar_angle != 0 && polar_angle != 1) {
      msg << "### FATAL ERROR in radiation class" << std::endl
          << "polar_angle must be 0 or 1" << std::endl;
      ATHENA_ERROR(msg);
    }
  }
}


// Assigns number of angles and octants for lab_davis and lab_magic
// Called both here and in src/mesh/meshblock.cpp
void NRRadiation::GetNumAnglesAndOctantsLab(int angle_flag, int ndim, int nmu,
    int &noct, int &nang) {
  if (angle_flag != ANGLE_FLAG_LAB_DAVIS && angle_flag != ANGLE_FLAG_LAB_MAGIC) {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetNumAnglesAndOctantsLab" << std::endl
        << "angle_flag must be " << ANGLE_FLAG_LAB_DAVIS << " or "
        << ANGLE_FLAG_LAB_MAGIC << ", got " << angle_flag << std::endl;
    ATHENA_ERROR(msg);
  }

  // Assign number of octants
  noct = std::pow(2, ndim);

  // Assign number of angles per octant
  int nang_per_oct = 0;

  if (ndim == 1) {
    nang_per_oct = nmu;

  } else if (ndim == 2) {
    if (angle_flag == ANGLE_FLAG_LAB_DAVIS) nang_per_oct = nmu * (nmu + 1) / 2;
    if (angle_flag == ANGLE_FLAG_LAB_MAGIC) nang_per_oct = nmu;

  } else if (ndim == 3) {
    if (angle_flag == ANGLE_FLAG_LAB_DAVIS) nang_per_oct = nmu * (nmu + 1) / 2;
    if (angle_flag == ANGLE_FLAG_LAB_MAGIC) nang_per_oct = nmu * nmu/2;
  }

  nang = nang_per_oct * noct;
}


// Assigns number of angles for rotating
// There are no octants: angle n = izeta*(2*npsi) + ipsi, with 2*nzeta zeta angles and
// 2*npsi psi angles (see AngularGrid(angle_flag, nzeta, npsi))
// Called both here and in src/mesh/meshblock.cpp
void NRRadiation::GetNumAnglesRotating(int angle_flag, int ndim, int polar_angle,
    int nzeta, int npsi, const char* coord, int &nang) {
  if (angle_flag != ANGLE_FLAG_COROTATING) {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetNumAnglesRotating" << std::endl
        << "angle_flag must be " << ANGLE_FLAG_COROTATING << ", got " << angle_flag
        << std::endl;
    ATHENA_ERROR(msg);
  }

  bool use_spherical_polar = (std::strcmp(coord, "spherical_polar") == 0);

  nang = 0;

  if (ndim == 1) {
    nang = 2*nzeta;

  } else if (ndim == 2) {
    if (use_spherical_polar) {
      nang = 4*nzeta*npsi;
    } else {
      // NOTE: AngularGrid only fills all of these angles for nzeta == 0 or npsi == 1
      if (npsi <= 1) nang = 4*nzeta; // only grid in theta
      if (nzeta == 0) nang = 2*npsi; // only grid in psi
      if (npsi > 1 && nzeta > 0) nang = 4*nzeta*npsi; // general case
    }

  } else if (ndim == 3) {
    nang = 4*nzeta*npsi;
  }

  // need to add two angles along the polar direction
  if (polar_angle)
    nang += 2;
}
