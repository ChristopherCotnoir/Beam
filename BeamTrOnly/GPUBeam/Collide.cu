
#include <cuComplex.h>
#include "Util.cu"


//!-------------------
//! STRONG-STRONG MODE
//!------------------------------------------------------------------

__device__
void
pRegimeStrongStrong(double &x, double &px, double &y, double &py,
		    double *dx_e, double *dxbar_e, double *dybar_e, double *dsig_x_e, double *dsig_y_e, int *dSe, double gamma_e, int Npart_e, double N_e,
		    double *dx_p, double *dxbar_p, double *dybar_p, double *dsig_x_p, double *dsig_y_p, int *dSp, double gamma_p, int Npart_p, double N_p,
		    double s, int curSlices, int sliceId){
  cuDoubleComplex eye = make_cuDoubleComplex(0.0 , 1.0);
  cuDoubleComplex z1, z2, w1, w2, fk;
  
  
  double xbar_e = dxbar_e[(curSlices - 1) - sliceId];
  double ybar_e = dybar_e[(curSlices - 1) - sliceId];

  double sig_x_e = dsig_x_e[(curSlices - 1) - sliceId];
  double sig_y_e = dsig_y_e[(curSlices - 1) - sliceId];
  int iSe = dSe[(curSlices - 1) - sliceId];

  double sig_x_p = dsig_x_p[sliceId];
  double sig_y_p = dsig_y_p[sliceId];
  int iSp = dSp[sliceId];


  double sxy_e = sig_x_e/sig_y_e;
  double syx_e = sig_y_e/sig_x_e;

  double sg_e = 1.0/sqrt(2.0 * (fabs(sig_x_e * sig_x_e - sig_y_e * sig_y_e)));
  double sg_p = 1.0/sqrt(2.0 * (fabs(sig_x_p * sig_x_p - sig_y_p * sig_y_p)));


  double x1 = x + s * px - xbar_e;
  double y1 = y + s * py - ybar_e;

  double eterm = exp(-0.50 * pow((x1/sig_x_e), 2.0) - 0.50 * pow((y1/sig_y_e), 2.0));
  double fcnst_e = (iSp*Re/gamma_e)*(2.0*sqrt(PI)*sg_p)*(N_p/Npart_p);
  double fcnst_p = (iSe*Rp/gamma_p)*(2.0*sqrt(PI)*sg_e)*(N_e/Npart_e);

  x1 = x1*sg_e;
  y1 = y1*sg_e;

  double x2 = syx_e*x1;
  double y2 = sxy_e*y1;
  double Fx = 0, Fy = 0;
  
  if (sig_x_e > sig_y_e){
    z1 = make_cuDoubleComplex(x1 + cuCreal(eye) * fabs(y1), cuCimag(eye)* fabs(y1));
    z2 = make_cuDoubleComplex(x2 + cuCreal(eye) * fabs(y2), cuCimag(eye)* fabs(y2));

    w1 = WOFZ(cuCreal(z1), cuCimag(z1));
    w2 = WOFZ(cuCreal(z2), cuCimag(z2));
        
    cuDoubleComplex tempC = cuCsub(w1, make_cuDoubleComplex(eterm*cuCreal(w2), eterm*cuCimag(w2)));
    fk = make_cuDoubleComplex(fcnst_p*cuCreal(tempC), fcnst_p*cuCimag(tempC));

    if (y1 > 0.0){
      Fx = cuCimag(fk);
      Fy = cuCreal(fk);
    }else{
      Fx = cuCimag(fk);
      Fy = -cuCreal(fk);
    }
  }else{
    z1 = make_cuDoubleComplex(y1 + cuCreal(eye) * fabs(x1), cuCimag(eye)* fabs(x1));
    z2 = make_cuDoubleComplex(y2 + cuCreal(eye) * fabs(x2), cuCimag(eye)* fabs(x2));

    w1 = WOFZ(cuCreal(z1), cuCimag(z1));
    w2 = WOFZ(cuCreal(z2), cuCimag(z2));
    
    cuDoubleComplex tempC = cuCsub(w1, make_cuDoubleComplex(eterm*cuCreal(w2), eterm*cuCimag(w2)));
    fk = make_cuDoubleComplex(fcnst_e*cuCreal(tempC), fcnst_e*cuCimag(tempC));
    if (x1 < 0.0){
      Fx = -cuCreal(fk);
      Fy = cuCimag(fk);
    }else{
      Fx = cuCreal(fk);
      Fy = cuCimag(fk);
    }    
  }
  
  
  x = x + s * Fx;
  px = px - Fx;
  y = y + s * Fy;
  py = py - Fy;
  

}

//!-------------------
//! STRONG-STRONG MODE
//!------------------------------------------------------------------
__device__
void
eRegimeStrongStrong(double &x, double &px, double &y, double &py,
		    double *dx_e, double *dxbar_e, double *dybar_e, double *dsig_x_e, double *dsig_y_e, int *dSe, double gamma_e, int Npart_e, double N_e,
		    double *dx_p, double *dxbar_p, double *dybar_p, double *dsig_x_p, double *dsig_y_p, int *dSp, double gamma_p, int Npart_p, double N_p,
		    double s, int curSlices, int sliceId){
  cuDoubleComplex eye = make_cuDoubleComplex(0.0 , 1.0);
  cuDoubleComplex z1, z2, w1, w2, fk;
  

  double sig_x_e = dsig_x_e[sliceId];
  double sig_y_e = dsig_y_e[sliceId];
  int iSe = dSe[sliceId];

  double xbar_p = dxbar_p[(curSlices - 1) - sliceId];
  double ybar_p = dybar_p[(curSlices - 1) - sliceId];
  double sig_x_p = dsig_x_p[(curSlices - 1) - sliceId];
  double sig_y_p = dsig_y_p[(curSlices - 1) - sliceId];
  int iSp = dSp[(curSlices - 1) - sliceId];

  double sxy_p = sig_x_p/sig_y_p;
  double syx_p = sig_y_p/sig_x_p;
  
  double sg_e = 1.0/sqrt(2.0 * (fabs(sig_x_e * sig_x_e - sig_y_e * sig_y_e)));
  double sg_p = 1.0/sqrt(2.0 * (fabs(sig_x_p * sig_x_p - sig_y_p * sig_y_p)));

  double fcnst_e = (iSp*Re/gamma_e)*(2.0*sqrt(PI)*sg_p)*(N_p/Npart_p);
  double fcnst_p = (iSe*Rp/gamma_p)*(2.0*sqrt(PI)*sg_e)*(N_e/Npart_e);

  double x1 = x - s * px - xbar_p;
  double y1 = y - s * py - ybar_p;
  double eterm = exp(-0.50*pow((x1/sig_x_p), 2.0) - 0.50 * pow((y1/sig_y_p), 2.0));

  x1 = x1*sg_p;
  y1 = y1*sg_p;

  double x2 = syx_p*x1;
  double y2 = sxy_p*y1;
  double Fx = 0, Fy = 0;

  
  
  if (sig_x_p > sig_y_p){
    z1 = make_cuDoubleComplex(x1 + cuCreal(eye) * fabs(y1), cuCimag(eye)* fabs(y1));
    z2 = make_cuDoubleComplex(x2 + cuCreal(eye) * fabs(y2), cuCimag(eye)* fabs(y2));
    
    w1 = WOFZ(cuCreal(z1), cuCimag(z1));
    w2 = WOFZ(cuCreal(z2), cuCimag(z2));
    
    
    cuDoubleComplex tempC = cuCsub(w1, make_cuDoubleComplex(eterm*cuCreal(w2), eterm*cuCimag(w2)));
    fk = make_cuDoubleComplex(fcnst_e*cuCreal(tempC), fcnst_e*cuCimag(tempC));

     if (y1 > 0.0){
       Fx = cuCimag(fk);
       Fy = cuCreal(fk);
     }else{
       Fx = cuCimag(fk);
       Fy = -cuCreal(fk);
     }
  }else{
    z1 = make_cuDoubleComplex(y1 + cuCreal(eye) * fabs(x1), cuCimag(eye)* fabs(x1));
    z2 = make_cuDoubleComplex(y2 + cuCreal(eye) * fabs(x2), cuCimag(eye)* fabs(x2));

    w1 = WOFZ(cuCreal(z1), cuCimag(z1));
    w2 = WOFZ(cuCreal(z2), cuCimag(z2));

    cuDoubleComplex tempC = cuCsub(w1, make_cuDoubleComplex(eterm*cuCreal(w2), eterm*cuCimag(w2)));
    fk = make_cuDoubleComplex(fcnst_e*cuCreal(tempC), fcnst_e*cuCimag(tempC));
    if (x1 < 0.0){
      Fx = -cuCreal(fk);
      Fy = cuCimag(fk);
    }else{
      Fx = cuCreal(fk);
      Fy = cuCimag(fk);
    }    
  }
  
  
  x = x - s * Fx;
  px = px - Fx;
  y = y - s * Fy;
  py = py - Fy;

}


__global__
void
applyKickGPU(double *dx_e, double *dxbar_e, double *dybar_e, double *dsig_x_e, double *dsig_y_e, int *dSe, double zmin_e, double gamma_e, int Npart_e, double N_e,
	     double *dx_p, double *dxbar_p, double *dybar_p, double *dsig_x_p, double *dsig_y_p, int *dSp, double zmin_p, double gamma_p, int Npart_p, double N_p,
	     double *dS, int numSlices, int curSlices, int iRegime, int isTopTriangle){
  int pid = blockIdx.x * blockDim.x + threadIdx.x;
  if(pid < Npart_e){
    double L_e    = 2.0 * fabs(zmin_e);  //! make the box size L = 2*|zmin|
    double dz_e  = L_e/(double)numSlices;
    int eSliceId = INT_MAX;

    double z_e  = dx_e[4 * Npart_e + pid];

    eSliceId  = (z_e - zmin_e)/dz_e + 1;
    eSliceId  = MIN(MAX(1,eSliceId), numSlices); //! equivalent to: if (iz>N) iz=N; if (iz<1) iz=1
    eSliceId--;
  
    double x  = dx_e[pid];
    double px = dx_e[Npart_e + pid];
    double y  = dx_e[2 * Npart_e + pid];
    double py = dx_e[3 * Npart_e + pid];    

    if(iRegime == 2 && eSliceId < curSlices && isTopTriangle){
      //printf("e %d\t%d\n", pid, eSliceId);
      double s = dS[eSliceId];
      eRegimeStrongStrong(x, px, y, py, 
			  dx_e, dxbar_e, dybar_e, dsig_x_e, dsig_y_e, dSe, gamma_e, Npart_e, N_e,
			  dx_p, dxbar_p, dybar_p, dsig_x_p, dsig_y_p, dSp, gamma_p, Npart_p, N_p,		    
			  s, curSlices, eSliceId);
    }

      
    if(iRegime == 2 && ((numSlices - 1 - eSliceId) < curSlices) && !isTopTriangle){
      double s = dS[eSliceId];
      eSliceId = eSliceId - (numSlices - curSlices);
      eRegimeStrongStrong(x, px, y, py, 
			  dx_e, dxbar_e, dybar_e, dsig_x_e, dsig_y_e, dSe, gamma_e, Npart_e, N_e,
			  dx_p, dxbar_p, dybar_p, dsig_x_p, dsig_y_p, dSp, gamma_p, Npart_p, N_p,		    
			  s, curSlices, eSliceId);
    }
    dx_e[pid] = x;
    dx_e[Npart_p + pid] = px;
    dx_e[2 * Npart_p + pid] = y;
    dx_e[3 * Npart_p + pid] = py;
  }

  if(pid < Npart_p){
    double L_p    = 2.0 * fabs(zmin_p);  //! make the box size L = 2*|zmin|
    double dz_p  = L_p/(double)numSlices;
    int pSliceId = INT_MAX;

    double z_p  = dx_p[4 * Npart_p + pid];

    pSliceId  = (z_p - zmin_p)/dz_p + 1;
    pSliceId  = MIN(MAX(1,pSliceId), numSlices); //! equivalent to: if (iz>N) iz=N; if (iz<1) iz=1
    pSliceId--;
  

    double x  = dx_p[pid];
    double px = dx_p[Npart_p + pid];
    double y  = dx_p[2 * Npart_p + pid];
    double py = dx_p[3 * Npart_p + pid];    

    //Top Triangle
    if(iRegime == 2 && pSliceId < curSlices && isTopTriangle){
      double s = dS[(curSlices - 1) - pSliceId];
      pRegimeStrongStrong(x, px, y, py, 
			  dx_e, dxbar_e, dybar_e, dsig_x_e, dsig_y_e, dSe, gamma_e, Npart_e, N_e,
			  dx_p, dxbar_p, dybar_p, dsig_x_p, dsig_y_p, dSp, gamma_p, Npart_p, N_p,		    
			  s, curSlices, pSliceId);
    }
    //Bottom Triangle
    if(iRegime == 2 && ((numSlices - 1 - pSliceId) < curSlices) && !isTopTriangle){
      double s = dS[(numSlices - 1) - pSliceId + (numSlices - curSlices)]; //Order of S is reverse of e-beam
      pSliceId = pSliceId - (numSlices - curSlices);
      pRegimeStrongStrong(x, px, y, py, 
			  dx_e, dxbar_e, dybar_e, dsig_x_e, dsig_y_e, dSe, gamma_e, Npart_e, N_e,
			  dx_p, dxbar_p, dybar_p, dsig_x_p, dsig_y_p, dSp, gamma_p, Npart_p, N_p,		    
			  s, curSlices, pSliceId);
      
    }
    dx_p[pid] = x;
    dx_p[Npart_p + pid] = px;
    dx_p[2 * Npart_p + pid] = y;
    dx_p[3 * Npart_p + pid] = py;  
  }

}


template<typename T>
__device__
void
Reduce(T *data, int N){
  // contiguous range pattern
  for(size_t offset = blockDim.x / 2; offset > 0; offset >>= 1){
    if(threadIdx.x < offset){
      for(int slice = 0; slice < N; ++slice){
	data[slice * blockDim.x + threadIdx.x] += data[slice * blockDim.x + threadIdx.x + offset];
      }
    }
    __syncthreads();
  }
  __syncthreads();
}


template<typename T>
__device__
void
init(T *data, int curSlices, T val){
  __syncthreads();
  
  for(int i = 0; i < curSlices; ++i){
    data[i * blockDim.x + threadIdx.x] = 0;
  }
}


//@param eORp Value of 1 implies e-beam and 2 implies p-beam  
__global__
void
sig_xy_reduce(double *dx, double *S, double *xbar, double *ybar, int Npart, 
	      double *opx, double zmin, int numSlices, int curSlices, int iRegime, int eORp, int isTopTriangle){
  extern __shared__ double sm[];
  int pid = blockIdx.x * blockDim.x + threadIdx.x;
  double  L    = 2.0 * fabs(zmin);  //! make the box size L = 2*|zmin|
  double dz    = L/(double)numSlices;
  int sliceId = INT_MAX;
  
  init(sm, curSlices, 0.0);

  double xval = 0, yval = 0;
  bool isParticleActive = false;
  if(pid < Npart){
    double x  = dx[pid];
    double px = dx[Npart + pid];
    double y  = dx[2 * Npart + pid];
    double py = dx[3 * Npart + pid];
    double z  = dx[4 * Npart + pid];

    sliceId  = (z - zmin)/dz + 1;
    sliceId      = MIN(MAX(1,sliceId), numSlices); //! equivalent to: if (iz>N) iz=N; if (iz<1) iz=1
    sliceId--;

    if(iRegime == 2 && sliceId < curSlices && eORp == 1 && isTopTriangle){
      double s = S[sliceId];
      xval = x - s*px - xbar[sliceId];
      yval = y - s*py - ybar[sliceId];
      isParticleActive = true;
    }

    if(iRegime == 2 && ((numSlices - 1 - sliceId) < curSlices) && eORp == 1 && !isTopTriangle){
      double s = S[sliceId];
      sliceId = sliceId - (numSlices - curSlices);
      xval = x - s*px - xbar[sliceId];
      yval = y - s*py - ybar[sliceId];
      isParticleActive = true;
    }
    
    if(iRegime == 2 && sliceId < curSlices && eORp == 2 && isTopTriangle){
      double s = S[(curSlices - 1) - sliceId];
      xval = x + s*px - xbar[sliceId];
      yval = y + s*py - ybar[sliceId];;
      isParticleActive = true;
    }
    if(iRegime == 2 && ((numSlices - 1 - sliceId) < curSlices) && eORp == 2 && !isTopTriangle){
      double s = S[(numSlices - 1) - sliceId + (numSlices - curSlices)]; //Order of S is reverse of e-beam
      sliceId = sliceId - (numSlices - curSlices);
      xval = x + s*px - xbar[sliceId];
      yval = y + s*py - ybar[sliceId];;
      isParticleActive = true;
      
    }
  }
  if(isParticleActive)
    sm[sliceId * blockDim.x + threadIdx.x] = xval * xval;

  __syncthreads();
  Reduce(sm, curSlices);
  if(threadIdx.x == 0){
    for(int i = 0; i < curSlices; ++i){
      opx[i * gridDim.x + blockIdx.x] = sm[i * blockDim.x];
    }
  }

  init(sm, curSlices, 0.0);
  if(isParticleActive)
    sm[sliceId * blockDim.x + threadIdx.x] = yval * yval;
  __syncthreads();
  Reduce(sm, curSlices);

  if(threadIdx.x == 0){
    for(int i = 0; i < curSlices; ++i){
      opx[(curSlices + i) * gridDim.x + blockIdx.x] = sm[i * blockDim.x];
    }
  }
}

__global__
void
xy_reduce_p(double *dx, double *S, double zmin,
	    int Npart, int numSlices, int curSlices, int iRegime,
	    double *opx, int *activeSlices, int isTopTriangle){
  extern __shared__ double sm[];
  int *intSm = (int *)&sm[0];
  int pid = blockIdx.x * blockDim.x + threadIdx.x;
  double  L    = 2.0 * fabs(zmin);  //! make the box size L = 2*|zmin|
  double dz    = L/(double)numSlices;
  int sliceId = INT_MAX;
  
  init(sm, curSlices, 0.0);

  double xval = 0, yval = 0;
  bool isParticleActive = false;
  if(pid < Npart){
    double x  = dx[pid];
    double px = dx[Npart + pid];
    double y  = dx[2 * Npart + pid];
    double py = dx[3 * Npart + pid];
    double z  = dx[4 * Npart + pid];

    sliceId  = (z - zmin)/dz + 1;
    sliceId      = MIN(MAX(1,sliceId), numSlices); //! equivalent to: if (iz>N) iz=N; if (iz<1) iz=1
    sliceId--;

    if(iRegime == 2 && sliceId < curSlices && isTopTriangle){
      double s = S[(curSlices - 1) - sliceId]; //Order of S is reverse of e-beam
      xval = x + s*px;
      yval = y + s*py;
      isParticleActive = true;
      //printf("%d %d %.16e\n", pid, sliceId, xval);
    }
    if(iRegime == 2 && ((numSlices - 1 - sliceId) < curSlices) && !isTopTriangle){
      double s = S[(numSlices - 1) - sliceId + (numSlices - curSlices)]; //Order of S is reverse of e-beam
      xval = x + s*px;
      yval = y + s*py;
      isParticleActive = true;
      //if(sliceId == 4)
      //printf("%d %d %.16e\n", pid, sliceId, xval);

      sliceId = sliceId - (numSlices - curSlices);
    }
  }
  if(isParticleActive)
    sm[sliceId * blockDim.x + threadIdx.x] = xval;

  __syncthreads();
  Reduce(sm, curSlices);
  if(threadIdx.x == 0){
    for(int i = 0; i < curSlices; ++i){
      opx[i * gridDim.x + blockIdx.x] = sm[i * blockDim.x];
    }
  }

  init(sm, curSlices, 0.0);
  if(isParticleActive)
    sm[sliceId * blockDim.x + threadIdx.x] = yval;
  __syncthreads();
  Reduce(sm, curSlices);

  if(threadIdx.x == 0){
    for(int i = 0; i < curSlices; ++i){
      opx[(curSlices + i) * gridDim.x + blockIdx.x] = sm[i * blockDim.x];
    }
  }

  //Count the number of particles in each slice
  init(intSm, curSlices, 0);
  if(isParticleActive){
    intSm[sliceId * blockDim.x + threadIdx.x] = 1;
  } 
  __syncthreads();
  Reduce(intSm, curSlices);
  
  if(threadIdx.x == 0){
    for(int i = 0; i < curSlices; ++i){
      activeSlices[i * gridDim.x + blockIdx.x] = intSm[i * blockDim.x];
    }
  }
}
__global__
void
xy_reduce_e(double *dx, double *S, double zmin,
	    int Npart, int numSlices, int curSlices, int iRegime,
	    double *opx, int *activeSlices, int isTopTriangle){
  extern __shared__ double sm[];
  int *intSm = (int *)&sm[0];
  int pid = blockIdx.x * blockDim.x + threadIdx.x;
  double  L    = 2.0 * fabs(zmin);  //! make the box size L = 2*|zmin|
  double dz    = L/(double)numSlices;
  int sliceId = INT_MAX;
  
  init(sm, curSlices, 0.0);

  double xval = 0, yval = 0;
  bool isParticleActive = false;
  if(pid < Npart){
    double x  = dx[pid];
    double px = dx[Npart + pid];
    double y  = dx[2 * Npart + pid];
    double py = dx[3 * Npart + pid];
    double z  = dx[4 * Npart + pid];

    sliceId  = (z - zmin)/dz + 1;
    sliceId      = MIN(MAX(1,sliceId), numSlices); //! equivalent to: if (iz>N) iz=N; if (iz<1) iz=1
    sliceId--;

    if(iRegime == 2 && sliceId < curSlices && isTopTriangle){   
      double s = S[sliceId];
      xval = x - s*px;
      yval = y - s*py;
      isParticleActive = true;
      //printf("%d %d %.16e\n", pid, sliceId, xval);
    }
    if(iRegime == 2 && ((numSlices - 1 - sliceId) < curSlices) && !isTopTriangle){
      double s = S[sliceId];
      sliceId = sliceId - (numSlices - curSlices);
      xval = x - s*px;
      yval = y - s*py;
      isParticleActive = true;
      //if(sliceId == 3)
      //printf("%d %d %.16e\n", pid, sliceId, xval);
    }
  }
  if(isParticleActive)
    sm[sliceId * blockDim.x + threadIdx.x] = xval;

  __syncthreads();
  Reduce<double>(sm, curSlices);
  if(threadIdx.x == 0){
    for(int i = 0; i < curSlices; ++i){
      opx[i * gridDim.x + blockIdx.x] = sm[i * blockDim.x];
    }
  }

  init(sm, curSlices, 0.0);
  if(isParticleActive)
    sm[sliceId * blockDim.x + threadIdx.x] = yval;
  __syncthreads();
  Reduce<double>(sm, curSlices);

  if(threadIdx.x == 0){
    for(int i = 0; i < curSlices; ++i){
      opx[(curSlices + i) * gridDim.x + blockIdx.x] = sm[i * blockDim.x];
    }
  }


  //Count the number of particles in each slice
  init(intSm, curSlices, 0);
  if(isParticleActive){
    //printf("%d %d\n", pid, sliceId);
    intSm[sliceId * blockDim.x + threadIdx.x] = 1;
  } 
  __syncthreads();
  Reduce(intSm, curSlices);
  
  if(threadIdx.x == 0){
    for(int i = 0; i < curSlices; ++i){
      activeSlices[i * gridDim.x + blockIdx.x] = intSm[i * blockDim.x];
    }
  }
}


//!----------------------------------------------------------------------------
//! Compute the mean and the standard deviation of the particle distribution 
//! in each of the two transversal coordinates for each beam 
//!----------------------------------------------------------------------------
void computeMeanAndSD(double *dx_e, double *dx_p, double *dS, 
		      int numKicks, double zmin_e, double zmin_p, 
		      double *&hxbar_e, double *&hybar_e, double *&hsig_x_e, double *&hsig_y_e, int *&hiSe,
		      double *&hxbar_p, double *&hybar_p, double *&hsig_x_p, double *&hsig_y_p, int *&hiSp,
		      BeamParams *bParams, int isTopTriangle){
  
  int numThreads = BLOCK_SIZE;
  int numBlocks = bParams->Npart_e/numThreads + ((bParams->Npart_e%numThreads)?1:0);

  double *dOpx = 0;
  int *dActiveSlices = 0;
  QuadDebug(cudaMalloc((void **)&dOpx, sizeof(double) * numBlocks * numKicks * 2));
  QuadDebug(cudaMalloc((void **)&dActiveSlices, sizeof(int) * numBlocks * numKicks));

  //Compute xbar_e, ybar_e by reducing x + S*px and y + S*py 

  //printf("\n\n");
  xy_reduce_e<<<numBlocks, numThreads, sizeof(double) * numKicks * numThreads>>>(dx_e, dS, zmin_e, bParams->Npart_e, bParams->N, numKicks, bParams->iRegime, dOpx, dActiveSlices, isTopTriangle);

  thrust::device_ptr<int> intPtr;
  thrust::device_ptr<double> dblPtr;
  
  std::cout << std::setprecision(18);
  std::cout << std::scientific;
  for(int i = 0; i < numKicks; ++i){
    intPtr = thrust::device_pointer_cast(dActiveSlices + i * numBlocks);
    hiSe[i] = thrust::reduce(intPtr, intPtr + numBlocks, (int) 0, thrust::plus<int>());
    //std::cout << hiSe[i] << "\n";

    dblPtr = thrust::device_pointer_cast(dOpx + i * numBlocks);    
    hxbar_e[i] = thrust::reduce(dblPtr, dblPtr + numBlocks, (double) 0.0, thrust::plus<double>());
    hxbar_e[i] /=(double)hiSe[i];
    //std::cout << hxbar_e[i] << "\n";

    dblPtr = thrust::device_pointer_cast(dOpx + numKicks * numBlocks + i * numBlocks);    
    hybar_e[i] = thrust::reduce(dblPtr, dblPtr + numBlocks, (double) 0, thrust::plus<double>())/(double)hiSe[i];
    //std::cout << hybar_e[i] << "\n";
  }

  double *dxbar_e = 0, *dybar_e = 0;
  QuadDebug(cudaMalloc((void **)&dxbar_e, sizeof(double) * numKicks));
  QuadDebug(cudaMalloc((void **)&dybar_e, sizeof(double) * numKicks));
  QuadDebug(cudaMemcpy(dxbar_e, hxbar_e, sizeof(double) * numKicks, cudaMemcpyHostToDevice));
  QuadDebug(cudaMemcpy(dybar_e, hybar_e, sizeof(double) * numKicks, cudaMemcpyHostToDevice));

  //Compute sig_xy by reducing
  sig_xy_reduce<<<numBlocks, numThreads, sizeof(double) * numKicks * numThreads>>>(dx_e, dS, dxbar_e, dybar_e, bParams->Npart_e, dOpx, zmin_e, bParams->N, numKicks, bParams->iRegime, 1, isTopTriangle);
  
  for(int i = 0; i < numKicks; ++i){
    dblPtr = thrust::device_pointer_cast(dOpx + i * numBlocks);    
    hsig_x_e[i] = sqrt(thrust::reduce(dblPtr, dblPtr + numBlocks, (double) 0, thrust::plus<double>())/(double)hiSe[i]);
    //std::cout << hsig_x_e[i] << "\n";

    dblPtr = thrust::device_pointer_cast(dOpx + numKicks * numBlocks + i * numBlocks);    
    hsig_y_e[i] = sqrt(thrust::reduce(dblPtr, dblPtr + numBlocks, (double) 0, thrust::plus<double>())/(double)hiSe[i]);
    //std::cout << hsig_y_e[i] << "\n";
  }

  
  numThreads = BLOCK_SIZE;
  numBlocks = bParams->Npart_p/numThreads + ((bParams->Npart_p%numThreads)?1:0);

  QuadDebug(cudaMalloc((void **)&dOpx, sizeof(double) * numBlocks * numKicks * 2));
  QuadDebug(cudaMalloc((void **)&dActiveSlices, sizeof(int) * numBlocks * numKicks));


  //For p-beam
  xy_reduce_p<<<numBlocks, numThreads, sizeof(double) * numKicks * numThreads>>>(dx_p, dS, zmin_p, bParams->Npart_p, bParams->N, numKicks, bParams->iRegime, dOpx, dActiveSlices, isTopTriangle);

  //display<double>(dOpx,1);
  
  for(int i = 0; i < numKicks; ++i){
    intPtr = thrust::device_pointer_cast(dActiveSlices + i * numBlocks);
    hiSp[i] = thrust::reduce(intPtr, intPtr + numBlocks, (int) 0, thrust::plus<int>());
    //std::cout << hiSp[i] << "\n";

    dblPtr = thrust::device_pointer_cast(dOpx + i * numBlocks);    
    hxbar_p[i] = thrust::reduce(dblPtr, dblPtr + numBlocks, (double) 0, thrust::plus<double>())/(double)hiSp[i];
    //std::cout << hxbar_p[i] << "\n";

    dblPtr = thrust::device_pointer_cast(dOpx + numKicks * numBlocks + i * numBlocks);    
    hybar_p[i] = thrust::reduce(dblPtr, dblPtr + numBlocks, (double) 0, thrust::plus<double>())/(double)hiSp[i];
    //std::cout << hybar_p[i] << "\n";

  }
  

  double *dxbar_p = 0, *dybar_p = 0;
  QuadDebug(cudaMalloc((void **)&dxbar_p, sizeof(double) * numKicks));
  QuadDebug(cudaMalloc((void **)&dybar_p, sizeof(double) * numKicks));
  QuadDebug(cudaMemcpy(dxbar_p, hxbar_p, sizeof(double) * numKicks, cudaMemcpyHostToDevice));
  QuadDebug(cudaMemcpy(dybar_p, hybar_p, sizeof(double) * numKicks, cudaMemcpyHostToDevice));

  //Compute sig_xy by reducing
  sig_xy_reduce<<<numBlocks, numThreads, sizeof(double) * numKicks * numThreads>>>(dx_p, dS, dxbar_p, dybar_p, bParams->Npart_p, dOpx, zmin_p, bParams->N, numKicks, bParams->iRegime, 2, isTopTriangle);

  //  display<double>(dOpx, 1);
  
  for(int i = 0; i < numKicks; ++i){
    dblPtr = thrust::device_pointer_cast(dOpx + i * numBlocks);    
    hsig_x_p[i] = sqrt(thrust::reduce(dblPtr, dblPtr + numBlocks, (double) 0, thrust::plus<double>())/(double)hiSp[i]);
    //std::cout << hsig_x_p[i] << "\n";

    dblPtr = thrust::device_pointer_cast(dOpx + numKicks * numBlocks + i * numBlocks);    
    hsig_y_p[i] = sqrt(thrust::reduce(dblPtr, dblPtr + numBlocks, (double) 0, thrust::plus<double>())/(double)hiSp[i]);
    //std::cout << hsig_y_p[i] << "\n";
  }
  

  QuadDebug(cudaFree(dActiveSlices));
  QuadDebug(cudaFree(dOpx));
}

void applyKick(double *dx_e, double *dx_p, double *hS, 
	       int numKicks, double zmin_e, double zmin_p, 
	       BeamParams *bParams, int isTopTriangle){

  double *xbar_e = new double[numKicks];
  double *ybar_e = new double[numKicks];
  int *iSe = new int[numKicks];
  double *sig_x_e = new double[numKicks];
  double *sig_y_e = new double[numKicks];

  double *xbar_p = new double[numKicks];
  double *ybar_p = new double[numKicks];
  int *iSp = new int[numKicks];
  double *sig_x_p = new double[numKicks];
  double *sig_y_p = new double[numKicks];

  double *dS = 0;
  QuadDebug(cudaMalloc((void **)&dS, sizeof(double) * bParams->N));
  QuadDebug(cudaMemcpy(dS, hS, sizeof(double) * bParams->N, cudaMemcpyHostToDevice));


  quad::timer::event_pair timer_node;
  quad::timer::start_timer(&timer_node);
  computeMeanAndSD(dx_e, dx_p, dS, numKicks, 
		   zmin_e, zmin_p, 
		   xbar_e, ybar_e, sig_x_e, sig_y_e, iSe, 
		   xbar_p, ybar_p, sig_x_p, sig_y_p, iSp, 
		   bParams, isTopTriangle);
  quad::timer::stop_timer(&timer_node, "ComputeMeanAndSD");
  
  std::cout << std::setprecision(18);
  std::cout << std::scientific;

  int nE = 0, nP = 0;
  for(int i = 0; i < numKicks; ++i){
    nE += iSe[i];
    nP += iSp[i];
    std::cout << "\n-------------------\n";
    std::cout << "slice\t" << i << "\n";
    std::cout << "eslice\t" << iSe[i] << "\n";
    std::cout << "pslice\t" << iSp[i] << "\n";
    std::cout << "xbar_e\t" << xbar_e[i] << "\n";
    std::cout << "ybar_e\t" << ybar_e[i] << "\n";
    std::cout << "xbar_p\t" << xbar_p[i] << "\n";
    std::cout << "ybar_p\t" << ybar_p[i] << "\n";
    std::cout << "sig_x_e\t" << sig_x_e[i] << "\n";
    std::cout << "sig_y_e\t" << sig_y_e[i] << "\n";
    std::cout << "sig_x_p\t" << sig_x_p[i] << "\n";
    std::cout << "sig_y_p\t" << sig_y_p[i] << "\n";
    
    //std::cout << hS[i] << "\n";
    //std::cout << xbar_e[i] << "\t" << ybar_e[i] << "\n";
    //std::cout << sig_x_e[i] << "\t" << sig_y_e[i] << "\n";
    //std::cout << xbar_p[i] << "\t" << ybar_p[i] << "\n";
    //std::cout << sig_x_p[i] << "\t" << sig_y_p[i] << "\n";
  }

  double *dxbar_e = 0, *dybar_e = 0, *dsig_x_e = 0, *dsig_y_e = 0;
  double *dxbar_p = 0, *dybar_p = 0, *dsig_x_p = 0, *dsig_y_p = 0;
  int *dSe = 0, *dSp = 0;

  QuadDebug(cudaMalloc((void **)&dxbar_e, sizeof(double) * numKicks));
  QuadDebug(cudaMalloc((void **)&dybar_e, sizeof(double) * numKicks));
  QuadDebug(cudaMalloc((void **)&dsig_x_e, sizeof(double) * numKicks));
  QuadDebug(cudaMalloc((void **)&dsig_y_e, sizeof(double) * numKicks));
  QuadDebug(cudaMalloc((void **)&dSe, sizeof(int) * numKicks));
  
  QuadDebug(cudaMalloc((void **)&dxbar_p, sizeof(double) * numKicks));
  QuadDebug(cudaMalloc((void **)&dybar_p, sizeof(double) * numKicks));
  QuadDebug(cudaMalloc((void **)&dsig_x_p, sizeof(double) * numKicks));
  QuadDebug(cudaMalloc((void **)&dsig_y_p, sizeof(double) * numKicks));
  QuadDebug(cudaMalloc((void **)&dSp, sizeof(int) * numKicks));

  QuadDebug(cudaMemcpy(dxbar_e, xbar_e, sizeof(double) * numKicks, cudaMemcpyHostToDevice));
  QuadDebug(cudaMemcpy(dybar_e, ybar_e, sizeof(double) * numKicks, cudaMemcpyHostToDevice));
  QuadDebug(cudaMemcpy(dsig_x_e, sig_x_e, sizeof(double) * numKicks, cudaMemcpyHostToDevice));
  QuadDebug(cudaMemcpy(dsig_y_e, sig_y_e, sizeof(double) * numKicks, cudaMemcpyHostToDevice));
  QuadDebug(cudaMemcpy(dSe, iSe, sizeof(int) * numKicks, cudaMemcpyHostToDevice));

  QuadDebug(cudaMemcpy(dxbar_p, xbar_p, sizeof(double) * numKicks, cudaMemcpyHostToDevice));
  QuadDebug(cudaMemcpy(dybar_p, ybar_p, sizeof(double) * numKicks, cudaMemcpyHostToDevice));
  QuadDebug(cudaMemcpy(dsig_x_p, sig_x_p, sizeof(double) * numKicks, cudaMemcpyHostToDevice));
  QuadDebug(cudaMemcpy(dsig_y_p, sig_y_p, sizeof(double) * numKicks, cudaMemcpyHostToDevice));
  QuadDebug(cudaMemcpy(dSp, iSp, sizeof(int) * numKicks, cudaMemcpyHostToDevice));

  int maxNpart = MAX(bParams->Npart_e, bParams->Npart_p);
  int numThreads = BLOCK_SIZE;
  int numBlocks = maxNpart/numThreads + ((maxNpart%numThreads)?1:0);
  
  quad::timer::start_timer(&timer_node);
  applyKickGPU<<<numBlocks, numThreads>>>(dx_e, dxbar_e, dybar_e, dsig_x_e, dsig_y_e, dSe, zmin_e, bParams->gamma_e, bParams->Npart_e, bParams->N_e,
					  dx_p, dxbar_p, dybar_p, dsig_x_p, dsig_y_p, dSp, zmin_p, bParams->gamma_p, bParams->Npart_p, bParams->N_p,
					  dS, bParams->N, numKicks, bParams->iRegime, isTopTriangle);
  quad::timer::stop_timer(&timer_node, "applyKickGPU");
  
  //display<double>(dx_p, bParams->Npart_e);
  //exit(1);
  QuadDebug(cudaFree(dS));
}

double slice(double *dx, int Npart){
  thrust::device_ptr<double> ptr;
  ptr = thrust::device_pointer_cast(dx);

  //! zmin (boundary) is the outlier
  double zmin = 9999;
  zmin = thrust::reduce(ptr + 4 * Npart, ptr + 5 * Npart, (double) zmin, thrust::minimum<double>());

  std::cout << std::setprecision(16);
  std::cout << "zmin\t:" << zmin << "\n";

  return zmin;   //! make the box size L = 2*|zmin|
}

void collide(double *dx_e, double *dx_p, BeamParams *bParams){
  double zmin_e = slice(dx_e, bParams->Npart_e);
  double zmin_p = slice(dx_p, bParams->Npart_p);
  double L_e = 2.0 * fabs(zmin_e);
  double L_p = 2.0 * fabs(zmin_p);

  //!-------------------------------------------------
  //! Set the drift length for multi-slice collision.
  //!-------------------------------------------------
  double delta_e = 0, delta_p = 0;
  int N = bParams->N;
  if (N == 1) {
    delta_e = 0.0;
    delta_p = 0.0;
  }else{
    delta_e = L_e/(double)N;
    delta_p = L_p/(double)N;
  }
  
  double *z_e = new double[N];
  double *z_p = new double[N];
  if (N == 1) {
    z_e[0] = 0.0;
    z_p[0] = 0.0;
  }else{
    int j = 0;
    for(int i = -(N - 1); i <= (N - 1); i+=2){
      j = j + 1;
      z_e[j - 1] = -i * delta_e/2.0;
      z_p[j - 1] = -i * delta_p/2.0;
    }
  }  


  //!-------------------------------------------------------
  //! Collide appropriate slices and apply drifts by drift
  //!-------------------------------------------------------
  //!   Example: N = 5
  //!--
  //! T  2:                                  11
  //! O  3:                          12              21
  //! P  4:                  13              22              31
  //!    5:          14              23              32              41
  //!__  6:  15              24              33              42              51
  //! B  7:          25              34              43              52
  //! O  8:                  35              44              53
  //! T  9:                          45              54
  //!   10:                                  55
  //!-----------------------------------------------------------------
  //! Top triangle
  //!---------------
  
  for(int i = 2; i <= N+1; ++i){
    double *S = new double[N];
    for(int j = 1; j <= i-1; ++j){
      S[j - 1] = (z_p[i-j-1]-z_e[j-1])/2.0;
    }
    quad::timer::event_pair timer_node;
    quad::timer::start_timer(&timer_node);
    
    applyKick(dx_e, dx_p, S, i - 1, zmin_e, zmin_p, bParams, 1);
    quad::timer::stop_timer(&timer_node, "applyKick");
  }
  
  //!-----------------
  //! Bottom triangle
  //!-----------------
  for(int i = N+2; i <= 2*N; ++i){
    double *S = new double[N];
    for(int j = i - N; j<= N; ++j){
      S[j - 1] = (z_p[i-j-1]-z_e[j - 1])/2.0;        
    }
    quad::timer::event_pair timer_node;
    quad::timer::start_timer(&timer_node);
    
    applyKick(dx_e, dx_p, S, 2*N - i + 1, zmin_e, zmin_p, bParams, 0);    
    quad::timer::stop_timer(&timer_node, "applyKick");
  }
  

}
