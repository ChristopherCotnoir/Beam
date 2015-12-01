
#define PI 3.14159265358979323846264338
#define Re 2.817940289458e-15    //! [meters]
#define Rp 1.534698e-18          //! [meters]
 
#define NCOL 6
#define PRECISION 2.2E-16

#define GCC_VERSION (__GNUC__ * 10000		      \
		     + __GNUC_MINOR__ * 100	      \
		     + __GNUC_PATCHLEVEL__)


struct BeamParams{
  int iTrackOnly;        //! flag for tracking only: 1: track <>1: 2-beam
  std::string gfEqns_file;//! gf equations file
  std::string Me_file;   //! matrix data file for the e-beam
  int Mef;               //! flag for matrix data file for the e-beam 
  int Norder_e;          //! order of the map for the e-beam
  int Npart_e;           //! number of simulation particles in the e-beam
  std::string ICe_file;  //! file containing ICs for the e-beam
  int ICef;              //! flag for file containing ICs for the e-beam 
  std::string Mp_file;   //! matrix data file for the e-beam
  int Mpf;               //! flag for matrix data file for the e-beam
  int Norder_p;          //! order of the map for the p-beam
  int Npart_p;           //! number of simulation particles in the p-beam
  std::string ICp_file;  //! file containing ICs for the p-beam
  int ICpf;              //! flag for file containing ICs for the p-beam 
  int Niter;            //! number of iterations
  int Nfreq;            //! frequency of dumping DF
  int NfreqLum;          //! frequency of dumping luminosity
  int iRegime;           //! type of beam-beam effect
  double sig_x0_e;       //! rms horizontal size of the e-beam
  double sig_y0_e;       //! rms vertical size of the e-beam
  double sig_z0_e;       //! rms size of the e-beam
  double sig_x0_p;       //! rms horizontal size of the p-beam
  double sig_y0_p;       //! rms vertical size of the p-beam
  double sig_z0_p;       //! rms size of the p-beam
  double sig_dE0;        //! rms energy spread
  double beta_x0;        //! horizontal beta* at IP
  double beta_y0;        //! vertical beta* at IP
  double N_e;               //! number of particles in the e-beam
  double N_p;               //! number of particles in the p-beam
  double E_e;            //! energy of the e-beam
  double E_p;            //! energy of the p-beam
  double f_0;            //! revolution frequency
  double nu_ex;          //! betatron tune in x for the e-beam
  double nu_ey;          //! betatron tune in y for the e-beam
  double nu_ez;          //! synchrotron tune for the e-beam
  double nu_px;          //! betatron tune in x for the p-beam
  double nu_py;          //! betatron tune in y for the p-beam
  double nu_pz;          //! synchrotron tune for the p-beam
  int N;                 //! number of slices in each beam
  int isGPU;		// flag for CPU/GPU execution: 1-GPU 0-CPU
  int isSymTr;		//0 - Tracking Only, 1 - Symplectic Tracking
  int icGen;            //! IC generation method
  int coord1;           //! First uniform lattice coordinate
  int coord2;           //! Second uniform lattice coordinate
  double x_l;            //! Uniform lattice lower bound for first coordinate
  double x_u;            //! Uniform lattice upper bound for first coordinate
  double y_l;            //! Uniform lattice lower bound for second coordinate
  double y_u;            //! Uniform lattice upper bound for second coordinate
  int coordNum1          //! Number of values for first uniform coordinate
  int coordNum2          //! Number of values for second uniform coordinate
  double gamma_e, gamma_p, Lc;
  double E_e0, E_p0;
  int jmaxord_e, jmaxord_p;

  int NSympFreq;
  bool strict_freq;
  bool log_in_background;
};

enum PARTICLE_TYPE {ELECTRON, PROTON};

