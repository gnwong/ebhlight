/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR FISHBONE-MONCRIEF TORUS                             *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#include "io.h"

// 
//  Notes:
//    * Use N3 == 1 at your own peril!
//    * If NO_BFIELD == 1, then we return an initial condition with B{123} = 0.
//    * If LOW_FLUX == 1, then we multiply q by   sin(M_PI/2. - th)  to flip the field across the midplane.
//    * I don't know why CLASSIC_BFIELD exists.
//
//  Basic structural flow:
//   (1) Generate FM torus in plasma. Requires rin/rmax be set.
//   (2) Sweep over fluid zones to find parameters needed for magnetic field initial condition.
//   (3) Set magnetic field A_phi and q per zone and generate magnetic field primitives.
//   (4) Update magnetic field primitives according to BHflux.
//
//  
//  The static int "MAD" should be thought of as a "magnetic-field-type" flag. You should not attempt to
//  read into the numbers. It's not worth it.
//   MAD == 0:  Penna+ (see MAD6) - like magnetic field initial condition a la BRR.
//   MAD == 1:  Standard SANE condition.
//   MAD == 2:  Standard MAD condition.
//   MAD == 3:  ??
//   MAD == 4:  BRR's newer semi-MAD condition where magnetic flux peaks around phi ~ 10-15 (/15, EHT convention).
//   MAD == 5:  ??
//   MAD == 6:  Penna, Kulkarni, Narayan 2013 AA 559 A116 (2013) magnetic field initial condition from A1.
//
//


#define CLASSIC_BFIELD (0)
#define LOW_FLUX (0)

// Use this to test the native Fishbone-Moncrief 
// torus. If (1), it will initialize the B^i:=0.
#define NO_BFIELD (0)

// Local functions
void coord_transform(double *Pr, int i, int j);
double lfish_calc(double rmax) ;

static int MAD;
static double BHflux;
static double beta;
static double u_jitter;
static double rin;
static double rmax;
static double Nloops; 


void set_problem_params()
{
  // register required parameters
  set_param("MAD", &MAD);
  set_param("BHflux", &BHflux);
  set_param("beta", &beta);
  
  // register optional parameters
  set_param_optional("u_jitter", &u_jitter);
  set_param_optional("rin", &rin);
  set_param_optional("rmax", &rmax);
  set_param_optional("Nloops", &Nloops);
 
  // set defaults for optional parameters
  u_jitter = 4.e-2;
  rin = 20.;
  rmax = 41.;
  Nloops = 8.;
}

void save_problem_params()
{
  hdf5_make_directory("macros");
  int tmp = CLASSIC_BFIELD;
  write_scalar(&tmp, "macros/CLASSIC_BFIELD", TYPE_INT);
  tmp = LOW_FLUX;
  write_scalar(&tmp, "macros/LOW_FLUX", TYPE_INT);
  tmp = NO_BFIELD;
  write_scalar(&tmp, "macros/NO_BFIELD", TYPE_INT);

  WRITE_HDR(MAD, TYPE_INT);
  WRITE_HDR(BHflux, TYPE_DBL);
  WRITE_HDR(beta, TYPE_DBL);
  WRITE_HDR(Nloops, TYPE_DBL);

  WRITE_HDR(u_jitter, TYPE_DBL);
  WRITE_HDR(rin, TYPE_DBL);
  WRITE_HDR(rmax, TYPE_DBL);
}

void init_prob()
{
  double r, th, sth, cth, ur, uh, up, u, rho, X[NDIM];
  struct of_geom *geom;

  // Disk interior
  double l, lnh, expm2chi, up1, DD, AA, SS, thin, sthin, cthin, DDin, AAin;
  double SSin, kappa, hm1;

  // Magnetic field
  static double A[N1+2*NG][N2+2*NG];
  double rho_av, rhomax, umax, bsq_ij, bsq_max, q;

  // Exists.
  if (N3 == 1) {
    fprintf(stderr, "! N3=1, defaulting to MAD=1 initial configuration.\n");
    MAD = 1; // SANE initial field not supported in 2D
    rin = 6.;
    rmax = 12.;
    Nloops = 1;
    DTd *= DTf;
    DTf = 1;
  }

  // Also exists.
  #if CLASSIC_BFIELD
  rin = 6.;
  rmax = 12.;
  #endif

  // ================================================================= //
  //                                                                   //
  // Generate Fishbone-Moncrief fluid state from rin, rmax conditions. //
  //                                                                   //
  // ================================================================= //

  l = lfish_calc(rmax);
  kappa = 1.e-3;

  rhomax = 0.;
  umax = 0.;
  ZSLOOP(-1, N1, -1, N2, -1, N3) {

    #ifdef RHOFL
    P[i][j][k][RHOFL] = 0.;
    P[i][j][k][UUFL] = 0.;
    #endif

    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    sth = sin(th);
    cth = cos(th);

    // Calculate lnh
    DD = r * r - 2. * r + a * a;
    AA = (r * r + a * a) * (r * r + a * a) -
        DD * a * a * sth * sth;
    SS = r * r + a * a * cth * cth;

    thin = M_PI / 2.;

    sthin = sin(thin);
    cthin = cos(thin);
    DDin = rin * rin - 2. * rin + a * a;
    AAin = (rin * rin + a * a) * (rin * rin + a * a)
        - DDin * a * a * sthin * sthin;
    SSin = rin * rin + a * a * cthin * cthin;

    if (r >= rin) {
      lnh =
          0.5 *
          log((1. +
         sqrt(1. +
              4. * (l * l * SS * SS) * DD / (AA *
                     sth *
                     AA *
                     sth)))
        / (SS * DD / AA))
          - 0.5 * sqrt(1. +
           4. * (l * l * SS * SS) * DD /
           (AA * AA * sth * sth))
          - 2. * a * r * l / AA -
          (0.5 *
           log((1. +
          sqrt(1. +
               4. * (l * l * SSin * SSin) * DDin /
               (AAin * AAin * sthin * sthin))) /
         (SSin * DDin / AAin))
           - 0.5 * sqrt(1. +
            4. * (l * l * SSin * SSin) *
            DDin / (AAin * AAin * sthin *
              sthin))
           - 2. * a * rin * l / AAin);
    } else {
      lnh = 1.;
    }

    // regions outside torus
    if (lnh < 0. || r < rin) {
      // Nominal values; real value set by fixup
      rho = 1.e-7 * RHOMIN;
      u = 1.e-7 * UUMIN;

      ur = 0.;
      uh = 0.;
      up = 0.;

      P[i][j][k][RHO] = rho;
      P[i][j][k][UU] = u;
      P[i][j][k][U1] = ur;
      P[i][j][k][U2] = uh;
      P[i][j][k][U3] = up;
    }
    /* region inside magnetized torus; u^i is calculated in
     * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
     * so it needs to be transformed at the end */
    else {
      hm1 = exp(lnh) - 1.;
      rho = pow(hm1 * (gam - 1.) / (kappa * gam),
          1. / (gam - 1.));
      u = kappa * pow(rho, gam) / (gam - 1.);
      ur = 0.;
      uh = 0.;

      // Calculate u^phi
      expm2chi = SS * SS * DD / (AA * AA * sth * sth);
      up1 =
          sqrt((-1. +
          sqrt(1. + 4. * l * l * expm2chi)) / 2.);
      up = 2. * a * r * sqrt(1. +
                 up1 * up1) / sqrt(AA * SS *
                 DD) +
          sqrt(SS / AA) * up1 / sth;


      P[i][j][k][RHO] = rho;
      if (rho > rhomax) rhomax = rho;
      P[i][j][k][UU] = u * (1. + u_jitter * (get_rand() - 0.5));
      if (u > umax && r > rin) umax = u;
      P[i][j][k][U1] = ur;
      P[i][j][k][U2] = uh;
      P[i][j][k][U3] = up;

      // Convert from 4-velocity to 3-velocity
      coord_transform(P[i][j][k], i, j);
    }

    P[i][j][k][B1] = 0.;
    P[i][j][k][B2] = 0.;
    P[i][j][k][B3] = 0.;
  } // ZSLOOP

  #if NO_BFIELD
  return;
  #endif

  // Normalize the densities so that max(rho) = 1. The actual normalization
  // is performed later so that we don't mess up the command flow of later
  // legacy code.
  umax = mpi_max(umax);
  rhomax = mpi_max(rhomax);


#if CLASSIC_BFIELD
  // Normalize densities
  ZSLOOP(-1, N1, -1, N2, -1, N3) {
    P[i][j][k][RHO] /= rhomax;
    P[i][j][k][UU] /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;
  fixup(P);
  bound_prim(P);
  //return;

  // Find vector potential at corners
  ZSLOOP(0, N1, 0, N2, 0, 0) A[i][j] = 0.;
  ZSLOOP(0, N1, 0, N2, 0, 0) {
    rho_av = 0.25*(P[i][j  ][0][RHO] + P[i-1][j  ][0][RHO] +
                   P[i][j-1][0][RHO] + P[i-1][j-1][0][RHO]);

    coord(i, j, k, CORN, X);
    bl_coord(X, &r, &th);

    q = rho_av/rhomax - 0.2;
    if (q > 0.) A[i][j] = q;
  }

  // Differentiate to find cell-centered B, and begin normalization
  bsq_max = 0.;
  ZLOOP {
    geom = get_geometry(i, j, k, CENT) ;

    // Flux-ct
    P[i][j][k][B1] = -(A[i][j] - A[i][j+1] + A[i+1][j] - A[i+1][j+1])/
                      (2.*dx[2]*geom->g);
    P[i][j][k][B2] = (A[i][j] + A[i][j+1] - A[i+1][j] - A[i+1][j+1])/
                     (2.*dx[1]*geom->g);
    P[i][j][k][B3] = 0.;

    bsq_ij = bsq_calc(P[i][j][k], geom);
    if (bsq_ij > bsq_max) bsq_max = bsq_ij;
  }
  bsq_max = mpi_max(bsq_max);

  // Normalize to set field strength
  double beta_act = (gam-1.)*umax/(0.5*bsq_max);
  double norm = sqrt(beta_act/beta);
  ZLOOP {
    P[i][j][k][B1] *= norm;
    P[i][j][k][B2] *= norm;
  }
#else

  // ============================================================================= //
  //                                                                               //                                                     
  // Now set up magentic field. This is a three-step process:                      //
  //                                                                               //
  //  (1) Run over fluid domain to find useful quantities like extent of the disk. //
  //  (2) Generate A_phi from q per zone.                                          //
  //  (3) Add extra BHflux at end.                                                 //
  //                                                                               //
  // ============================================================================= //

  // Penna+ 2013 - like magnetic field
  if (MAD==0 || MAD==6) {

    // these numbers are used to set the extent of the magnetic 
    // field within the disk. magnetic field is set to be non-
    // zero within the region (in the midplane) where density
    // rho > BRHOMIN. old values were 23->374 for 20,41 large 
    // disk and 12->84 for 10,21 medium disk.
    double BRHOMIN_inside = 2.e-1;
    double BRHOMIN_outside = 2.e-2;
    double rstart = 1.e10;
    double rend = 0.;

    // find rstart and rend. normalize along the way.
    ZSLOOP(-1, N1, -1, N2, -1, N3) {

      // normalize rho and uu
      P[i][j][k][RHO] /= rhomax;
      P[i][j][k][UU] /= rhomax;

      // identify rstart and rend by rho > 0.2 condition. only
      // run this on true midplane zones.
      int jglobal = j - NG + global_start[2];
      if (jglobal == N2TOT/2) {
        if (rstart==1.e10) {
          if (P[i][j][k][RHO]>BRHOMIN_inside) {
            coord(i,j,k, FACE1, X);
            bl_coord(X, &rstart, &th);
          }
        } else if (rend==0. && P[i][j][k][RHO]<BRHOMIN_outside) {
          coord(i,j,k, CENT, X);
          bl_coord(X, &rend, &th);
        }
      }
    }

    rstart = mpi_min(rstart);
    rend = mpi_max(rend);

    // print some diagnostic information ot the terminal
    if (mpi_io_proc()) {
      fprintf(stderr, "Problem file: torus\n  MAD: %d\n", MAD);
      fprintf(stderr, "  Nloops: %d\n  rstart: %g\n", (int)Nloops, rstart);
      fprintf(stderr, "  rend: %g\n", rend);
      fprintf(stderr, "\n");
    }
    
    // update rho, uu max values now that they're normalized.
    umax /= rhomax;
    rhomax = 1.;
    fixup(P);
    bound_prim(P);

    // now that we have rstart and rend, and that rho & uu are
    // properly normalized, find u_cm and u_m[...]
    double u_cm = 0.;
    double u_m[N1TOT];
    for (int i=0; i<N1TOT; ++i) u_m[i] = 0.;

    ZSLOOP(0, N1, 0, N2, 0, 0) {
      
      int iglobal = i - NG + global_start[1];
      int jglobal = j - NG + global_start[2];
      if (jglobal == N2TOT/2) {
        coord(i,j,k, FACE1, X);
        bl_coord(X, &r, &th);
        if ( r > rend && u_cm == 0. ) {
          u_cm = P[i][j][k][UU];
        }
        u_m[iglobal] = P[i][j][k][UU];
      }

    }

    // synchronize u_cm and u_m over all ranks.
    u_cm = mpi_max(u_cm);
    for (int i=0; i<N1TOT; ++i) {
      u_m[i] = mpi_max(u_m[i]);
    }

    // first find corner-centered vector potential
    ZSLOOP(0, N1, 0, N2, 0, 0) A[i][j] = 0.;
    ZSLOOP(0, N1, 0, N2, 0, 0) {

      rho_av = 0.25*(P[i][j  ][0][RHO] + P[i-1][j  ][0][RHO] +
                     P[i][j-1][0][RHO] + P[i-1][j-1][0][RHO])
        *(1. + 0.0*(get_rand() - 0.5));

      double uu_av = 0.25*(P[i][j  ][0][UU] + P[i-1][j  ][0][UU] +
                           P[i][j-1][0][UU] + P[i-1][j-1][0][UU]);

      coord(i, j, k, CORN, X);
      bl_coord(X, &r, &th);

      q = 0.;
      if (MAD == 0) {
        // this is the BRR alternative
        q = r/rstart*rho_av/rhomax - 0.2;
      } else if (MAD == 6) {
        int iglobal = i - NG + global_start[1];
        double Ucm = u_m[iglobal] - u_cm;
        double Uc = uu_av - u_cm;
        q = (Uc/Ucm - 0.2) / 0.8;
        q *= pow(sin(th), 3.);
      }

      double lamb = (log(rend) - log(rstart))/(M_PI*Nloops);
      A[i][j] = 0.;
      if (r > rstart && r < rend && q > 0.) {
        if (N3 > 1) {
          A[i][j] = q*sin((log(r) - log(rstart))/lamb);
          #if LOW_FLUX
          A[i][j] *= sin(M_PI/2. - th);
          #endif
        } else {
          A[i][j] = q;
        }
      }
    }

    // differentiate to find cell-centered B, compute bsq_max and beta_max along
    // the way.
    bsq_max = 0.;
    ZLOOP {
      geom = get_geometry(i, j, k, CENT) ;

      // Flux-ct
      P[i][j][k][B1] = -(A[i][j] - A[i][j + 1]
          + A[i + 1][j] - A[i + 1][j + 1]) /
          (2. * dx[2] * geom->g);
      P[i][j][k][B2] = (A[i][j] + A[i][j + 1]
               - A[i + 1][j] - A[i + 1][j + 1]) /
               (2. * dx[1] * geom->g);

      P[i][j][k][B3] = 0.;

      bsq_ij = bsq_calc(P[i][j][k], geom);
      if (bsq_ij > bsq_max) bsq_max = bsq_ij;
    }
    bsq_max = mpi_max(bsq_max);

    // ================================================================ //
    // Normalization happens here. One method is to bound beta, another //
    // is to try to force the total flux to be equal per loop.          //
    // ===============================================================  //

    double ldrloop = (log(rend) - log(rstart))/Nloops;

    if (MAD == 0) {
      // bound by plasma beta
      double facs[(int)Nloops];

      for (int n = 0; n < Nloops; n++) {
        double rmin = exp(log(rstart) + n*ldrloop);
        double rmax = exp(log(rstart) + (n+1)*ldrloop);
        double beta_min = 1.e100;
        double bsq_max = 0.;
        double umax = 0.;
        ZLOOP {
          coord(i, j, k, CENT, X);
          bl_coord(X, &r, &th);
          if (r > rmin && r < rmax) {
            // Inside loop
            geom = get_geometry(i, j, k, CENT);
            double bsq = bsq_calc(P[i][j][k], geom);
            if (bsq > bsq_max) bsq_max = bsq;
            if (P[i][j][k][UU] > umax) umax = P[i][j][k][UU];
            double beta_ij = (gam - 1.)*P[i][j][k][UU]/(0.5*bsq);
            if (beta_ij < beta_min) {
              beta_min = beta_ij;
            }
          }
        }

        if (mpi_nprocs() == 1) {
          facs[n] = sqrt(beta_min/beta);
          printf("fac[%i] = %e\n", n, sqrt(beta_min/beta));
        }

        // Rescale vector potential for this loop
        ZSLOOP(0,N1,0,N2,0,0) {
          coord(i, j, k, CORN, X);
          bl_coord(X, &r, &th);
          if (r > rmin && r < rmax) {
            if (mpi_nprocs() > 1) {
              A[i][j] *= facs[n];
            } else {
              A[i][j] *= sqrt(beta_min/beta);
            }
          }
        }
      }

    } else if (MAD == 6) {

      double maxPhi = 0.;
      double init_fluxes[(int)Nloops];
      double min_betas[(int)Nloops];

      // try to keep fluxes the same in each loop
      for (int n=0; n<Nloops; ++n) {

        // set beta absurdly high
        min_betas[n] = 1.e100;

        // get loop geometry
        double rmin = exp(log(rstart) + n*ldrloop);
        double rmax = exp(log(rstart) + (n+1)*ldrloop);

        // calculate total flux inside of this loop and get min beta
        double Phi_scaled = 0.;
        ISLOOP(5, N1-1) {
          JSLOOP(0, N2-1) {
            int jglobal = j - NG + global_start[2];
            int k = NG;

            // only compute for midplane
            if (jglobal == N2TOT/2) {

              coord(i, j, k, CENT, X);
              bl_coord(X, &r, &th);

              // only consider in the loop
              if (r > rmin && r < rmax) {

                // compute flux through midplane (field symmetry forces orthogonality)
                double B2net =  (A[i][j] + A[i][j+1] - A[i+1][j] - A[i+1][j+1])/
                                (2.*dx[1]*ggeom[i][j][CENT].g);
                Phi_scaled += ggeom[i][j][CENT].g*2.*M_PI*dx[1]*fabs(B2net);
                
                // compute beta (we set P before this normalization step)
                geom = get_geometry(i, j, k, CENT) ;
                double bsq_ij = bsq_calc(P[i][j][k], geom);
                double beta_ij = (gam - 1.)*P[i][j][k][UU]/(0.5*bsq_ij);
                if (beta_ij < min_betas[n]) min_betas[n] = beta_ij;
              }
            }
          }
        }

        init_fluxes[n] = mpi_reduce(Phi_scaled);
        min_betas[n] = mpi_min(min_betas[n]);
        if (init_fluxes[n] > maxPhi) maxPhi = init_fluxes[n];
      }

      // now actually normalize
      for (int n=0; n<Nloops; ++n) {

        // a kludge, to be sure, but a welcome(?) one.
        double sizefact = 20.;
        if (rin < 15) sizefact = 4.;

        double rmin = exp(log(rstart) + n*ldrloop);
        double rmax = exp(log(rstart) + (n+1)*ldrloop);
        double normalization = maxPhi / init_fluxes[n] * sizefact;
        // the extra factor at the end of the normalization can be
        // used to adjust the overall (e.g.) beta everywhere. this
        // is nomially 1. for medium disks and 20. for large disks
        // (10,20 vs 20,41). decreasing this factor will result in 
        // an increase in the overall minimum beta per loop.
        
        ZSLOOP(0,N1,0,N2,0,0) {
          coord(i, j, k, CENT, X);
          bl_coord(X, &r, &th);
          if (r > rmin && r < rmax) {
            A[i][j] *= normalization; 
          }
        }
          
      }
    }

    // ==================================== //
    // Housekeeping. Clean up B primitives. //
    // ==================================== //
    ZLOOP {
      geom = get_geometry(i, j, k, CENT) ;
      // Flux-ct
      P[i][j][k][B1] = -(A[i][j] - A[i][j + 1]
          + A[i + 1][j] - A[i + 1][j + 1]) /
          (2. * dx[2] * geom->g);
      P[i][j][k][B2] = (A[i][j] + A[i][j + 1]
               - A[i + 1][j] - A[i + 1][j + 1]) /
               (2. * dx[1] * geom->g);
    }
  } 


  if (MAD != 0 && MAD != 6) {  // begin MAD = {1,2,3,4,5} block

    // normalize density and internal energy. this is done here so as not to mess with 
    // previous definitions and code flow
    ZSLOOP(0, N1, 0, N2, 0, 0) {
      P[i][j][k][RHO] /= rhomax;
      P[i][j][k][UU] /= rhomax;
    }
    umax /= rhomax;
    rhomax = 1.;
    fixup(P);
    bound_prim(P);

    // now we start to set up the magnetic field by first finding the corner-centered
    // vector potential.
    ZSLOOP(0, N1, 0, N2, 0, 0) A[i][j] = 0.;
    ZSLOOP(0, N1, 0, N2, 0, 0) {
      rho_av = 0.25*(P[i][j  ][0][RHO] + P[i-1][j  ][0][RHO] +
                     P[i][j-1][0][RHO] + P[i-1][j-1][0][RHO])
        *(1. + 0.0*(get_rand() - 0.5));

      coord(i, j, k, CORN, X);
      bl_coord(X, &r, &th);

      if (N3 > 1) {

        if (MAD == 1) {
          q = rho_av/rhomax - 0.2;
        } else if (MAD == 2) {
          q = pow(sin(th),3)*pow(r/rin,3.)*exp(-r/400)*rho_av/rhomax - 0.2;
        } else if (MAD == 3) {
          q = (r/rin)*rho_av/rhomax - 0.2;
        } else if (MAD == 4) {
          q = pow(r/rin,2)*rho_av/rhomax - 0.2;
        } else if (MAD == 5) {
          q = pow(r/rin,3)*rho_av/rhomax - 0.2;
        } else {
          printf("MAD = %i not supported!\n", MAD);
          exit(-1);
        }

      } else {
        q = rho_av/rhomax - 0.2;
      }

      A[i][j] = 0.;
      if (q > 0.) A[i][j] = q;
    }

    // given the vector potential A_phi, differentiate to find the cell-centered
    // values for B. along the way, compute bsq_max for normalization later on.
    // recall that this normalization is performed to restrict plasma beta into 
    // a certain range based on the runtime parameter file.
    bsq_max = 0.;
    ZLOOP {
      geom = get_geometry(i, j, k, CENT);

      // Flux-ct
      P[i][j][k][B1] = -(A[i][j] - A[i][j + 1]
          + A[i + 1][j] - A[i + 1][j + 1]) /
          (2. * dx[2] * geom->g);
      P[i][j][k][B2] = (A[i][j] + A[i][j + 1]
               - A[i + 1][j] - A[i + 1][j + 1]) /
               (2. * dx[1] * geom->g);

      P[i][j][k][B3] = 0.;

      bsq_ij = bsq_calc(P[i][j][k], geom);
      if (bsq_ij > bsq_max) bsq_max = bsq_ij;
    }
    bsq_max = mpi_max(bsq_max);

    // normalize to set field strength according to input beta parameter.
    double beta_act = (gam - 1.)*umax/(0.5*bsq_max);
    double norm = sqrt(beta_act/beta);
    ZLOOP {
      P[i][j][k][B1] *= norm;
      P[i][j][k][B2] *= norm;
    }
  }  // end MAD = {1,2,3,4,5} block


  // ================================================================= //
  //                                                                   //
  // Use this section to add an additional net magnetic field into the //
  // torus. The strength of this field is set by the BHflux runtime    //
  // parameter. If BHflux == 0, you can ignore this code.              //
  //                                                                   //
  // ================================================================= //

  // Initialize a net magnetic field inside the initial torus
  ZSLOOP(0, N1, 0, N2, 0, 0) {
    A[i][j] = 0.;
    coord(i,j,k,CORN,X);
    bl_coord(X, &r, &th);
    double x = r*sin(th);
    double z = r*cos(th);
    double a_hyp = rin;
    double b_hyp = 3.*rin;
    double x_hyp = a_hyp*sqrt(1. + pow(z/b_hyp,2));

    q = (pow(x,2) - pow(x_hyp,2))/pow(x_hyp,2);
    if (x < x_hyp) {
      A[i][j] = 10.*q;
    }
  }

  // Evaluate net flux
  double Phi = 0.;
  ISLOOP(5, N1-1) {
    JSLOOP(0, N2-1) {
      int jglobal = j - NG + global_start[2];
      int k = NG;
      if (jglobal == N2TOT/2) {
        coord(i, j, k, CENT, X);
        bl_coord(X, &r, &th);
        if (r < rin) {
          double B2net =  (A[i][j] + A[i][j+1] - A[i+1][j] - A[i+1][j+1])/
                          (2.*dx[1]*ggeom[i][j][CENT].g);
          Phi += ggeom[i][j][CENT].g*2.*M_PI*dx[1]*fabs(B2net)/N3CPU;
        }
      }
    }
  }

  if (global_start[1] == 0) {
    JSLOOP(0, N2/2-1) {
      int i = 5 + NG;
      int k = NG;
      coord(i, j, k, CENT, X);
      bl_coord(X, &r, &th);
      double B1net = -(A[i][j] - A[i][j+1] + A[i+1][j] - A[i+1][j+1])/(2.*dx[2]*ggeom[i][j][CENT].g);
      Phi += ggeom[i][j][CENT].g*dx[2]*2.*M_PI*fabs(B1net)/N3CPU;
    }
  }
  Phi = mpi_reduce(Phi);

  double norm = BHflux/(Phi + SMALL);

  ZLOOP {
    geom = get_geometry(i, j, k, CENT) ;

    // Flux-ct
    P[i][j][k][B1] += -norm*(A[i][j] - A[i][j + 1]
        + A[i + 1][j] - A[i + 1][j + 1]) /
        (2. * dx[2] * geom->g);
    P[i][j][k][B2] += norm*(A[i][j] + A[i][j + 1]
             - A[i + 1][j] - A[i + 1][j + 1]) /
             (2. * dx[1] * geom->g);
  }

#endif // CLASSIC_BFIELD
}

// Convert Boyer-Lindquist four-velocity to MKS 3-velocity
void coord_transform(double *Pr, int ii, int jj)
{
  double X[NDIM], r, th, ucon[NDIM], trans[NDIM][NDIM], tmp[NDIM];
  double AA, BB, CC, discr;
  double alpha, gamma, beta[NDIM];
  struct of_geom *geom, blgeom;

  coord(ii, jj, 0, CENT, X);
  bl_coord(X, &r, &th);
  blgset(ii, jj, &blgeom);

  ucon[1] = Pr[U1];
  ucon[2] = Pr[U2];
  ucon[3] = Pr[U3];

  AA = blgeom.gcov[0][0];
  BB = 2. * (blgeom.gcov[0][1] * ucon[1] +
       blgeom.gcov[0][2] * ucon[2] +
       blgeom.gcov[0][3] * ucon[3]);
  CC = 1. +
      blgeom.gcov[1][1] * ucon[1] * ucon[1] +
      blgeom.gcov[2][2] * ucon[2] * ucon[2] +
      blgeom.gcov[3][3] * ucon[3] * ucon[3] +
      2. * (blgeom.gcov[1][2] * ucon[1] * ucon[2] +
      blgeom.gcov[1][3] * ucon[1] * ucon[3] +
      blgeom.gcov[2][3] * ucon[2] * ucon[3]);

  discr = BB * BB - 4. * AA * CC;
  ucon[0] = (-BB - sqrt(discr)) / (2. * AA);
  // This is ucon in BL coords

  // transform to Kerr-Schild
  // Make transform matrix
  memset(trans, 0, 16*sizeof(double));
  for (int mu = 0; mu < NDIM; mu++) {
    trans[mu][mu] = 1.;
  }
  trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
  trans[3][1] = a / (r * r - 2. * r + a * a);

  // Transform from BL to KS coordinates
  for (int mu = 0; mu < NDIM; mu++) tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
     tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) ucon[mu] = tmp[mu];

  // Transform from KS to MKS coordinates
  double invtrans[NDIM][NDIM];
  set_dxdX(X, invtrans);
  invert(&invtrans[0][0], &trans[0][0]);

  for (int mu = 0; mu < NDIM; mu++) tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
     tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) ucon[mu] = tmp[mu];

  //ucon[1] *= (1. / (r - R0));
  //ucon[1] /= dr_dx(X[1]);
  //ucon[2] *=
  //    (1. / (M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2])));
  //ucon[2] /= dth_dx(X[2]);

  // Solve for v. Use same u^t, unchanged under KS -> KS'
  geom = get_geometry(ii, jj, 0, CENT) ;
  alpha = 1.0 / sqrt( -geom->gcon[0][0] ) ;
  gamma = ucon[0] * alpha;

  beta[1] = alpha * alpha * geom->gcon[0][1];
  beta[2] = alpha * alpha * geom->gcon[0][2];
  beta[3] = alpha * alpha * geom->gcon[0][3];

  Pr[U1] = ucon[1] + beta[1] * gamma / alpha;
  Pr[U2] = ucon[2] + beta[2] * gamma / alpha;
  Pr[U3] = ucon[3] + beta[3] * gamma / alpha;
}

double lfish_calc(double r)
{
  return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
     ((-2. * a * r *
       (pow(a, 2) - 2. * a * sqrt(r) +
        pow(r,
      2))) / sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
      ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) + pow(a, 2) *
      (2. + r))) / sqrt(1 + (2. * a) / pow (r, 1.5) - 3. / r)))
    / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
       (pow(a, 2) + (-2. + r) * r))
      );
}

