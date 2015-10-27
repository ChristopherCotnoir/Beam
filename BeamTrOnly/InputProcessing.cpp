#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include "common.h"
#include "Beam.h"
#include "IO.h"


IO::IO(){
  inputConfig["input"] = "input/params.sim";
  paths["input"] = "input/";
  paths["output"] = "output/";
  
  outputConfig["IC_e"] = "output/dump.IC_e";
  outputConfig["dump.ebeam"] = "output/dump.ebeam";

  outputConfig["IC_p"] = "output/dump.IC_p";
  outputConfig["dump.pbeam"] = "output/dump.pbeam";

  outputConfig["luminosity"] = "output/luminosity.out";
    
  for (std::map<std::string, std::string>::iterator it = outputConfig.begin(); it != outputConfig.end(); ++it){
    std::ofstream file(it->second.c_str(), std::ios::out | std::ios::trunc);
  }
}


InputProcessing::InputProcessing(int pid){
  this->pending_log = false;
  this->log_in_progress = false;
  this->pid = pid;
}

void 
InputProcessing::ReadInput(BeamParams* bParams){
  std::ifstream beamInput(inputConfig["input"].c_str());
  if (!beamInput.is_open()){
    std::string msg = "Cannot open file " + inputConfig["input"];
    Abort(msg.c_str());
  }
  std::string line;

  beamInput >> bParams->iTrackOnly;getline(beamInput, line);
  beamInput >> bParams->Me_file;getline(beamInput, line);
  beamInput >> bParams->Mef;getline(beamInput, line);
  beamInput >> bParams->Norder_e;getline(beamInput, line);
  beamInput >> bParams->Npart_e;getline(beamInput, line);
  beamInput >> bParams->ICe_file;getline(beamInput, line);
  beamInput >> bParams->ICef;getline(beamInput, line);
  beamInput >> bParams->Mp_file;getline(beamInput, line);
  beamInput >> bParams->Mpf;getline(beamInput, line);
  beamInput >> bParams->Norder_p;getline(beamInput, line);
  beamInput >> bParams->Npart_p;getline(beamInput, line);
  beamInput >> bParams->ICp_file;getline(beamInput, line);
  beamInput >> bParams->ICpf;getline(beamInput, line);
  beamInput >> bParams->Niter;getline(beamInput, line);
  beamInput >> bParams->Nfreq;getline(beamInput, line);
  beamInput >> bParams->NfreqLum;getline(beamInput, line);
  beamInput >> bParams->iRegime;getline(beamInput, line);
  beamInput >> bParams->sig_x0_e;getline(beamInput, line);
  beamInput >> bParams->sig_y0_e;getline(beamInput, line);
  beamInput >> bParams->sig_z0_e;getline(beamInput, line);
  beamInput >> bParams->sig_x0_p;getline(beamInput, line);
  beamInput >> bParams->sig_y0_p;getline(beamInput, line);
  beamInput >> bParams->sig_z0_p;getline(beamInput, line);
  beamInput >> bParams->sig_dE0;getline(beamInput, line);
  beamInput >> bParams->beta_x0;getline(beamInput, line);
  beamInput >> bParams->beta_y0;getline(beamInput, line);
  beamInput >> bParams->N_e;getline(beamInput, line);
  beamInput >> bParams->N_p;getline(beamInput, line);
  beamInput >> bParams->E_e;getline(beamInput, line);
  beamInput >> bParams->E_p;getline(beamInput, line);
  beamInput >> bParams->f_0;getline(beamInput, line);
  beamInput >> bParams->nu_ex;getline(beamInput, line);
  beamInput >> bParams->nu_ey;getline(beamInput, line);
  beamInput >> bParams->nu_ez;getline(beamInput, line);
  beamInput >> bParams->nu_px;getline(beamInput, line);
  beamInput >> bParams->nu_py;getline(beamInput, line);
  beamInput >> bParams->nu_pz;getline(beamInput, line);
  beamInput >> bParams->N;getline(beamInput, line);

  beamInput >> bParams->gfEqns_file;getline(beamInput, line);

  beamInput >> bParams->x_bound;getline(beamInput, line);
  beamInput >> bParams->y_bound;getline(beamInput, line);
  beamInput >> bParams->NSympFreq;getline(beamInput, line);
  
  beamInput >> bParams->log_in_background;getline(beamInput, line);
  beamInput >> bParams->strict_freq;getline(beamInput, line);
  
  beamInput >> bParams->isGPU;getline(beamInput, line);
  beamInput >> bParams->isSymTr;getline(beamInput, line);
 
  //!------------------------------------------------
  //! Compute the parameters used in the simulation.
  //!------------------------------------------------
  bParams->E_e0 = 511e3;                //! [eV]
  bParams->E_p0 = 938e6;                //! [eV]
  bParams->gamma_e = bParams->E_e/bParams->E_e0;
  bParams->gamma_p = bParams->E_p/bParams->E_p0;
  bParams->Lc      = 1e-04*bParams->N_e*bParams->N_p*bParams->f_0/(2.0*PI);         //! [cm^{-2} s^{-1}]

  beamInput.close();
}

void 
InputProcessing::ReadMap(std::string mapFileName, int *Nrow, double *M, int *it, int Norder, int &jmaxOrd){
  std::string mapFilePath = paths["input"] + mapFileName;
  std::ifstream mapFile(mapFilePath.c_str());
  if (!mapFile.is_open()){
    std::string msg = "Cannot open file " + mapFileName;
    Abort(msg.c_str());
  }
  std::string line;
  int iTmp = 0, iOrder, iTmp_prev = 0;
  double Mtmp1;
  int ittmp[6];
  int j = 0, i = 1;
  while(!mapFile.eof()){
    mapFile >> iTmp >> Mtmp1 >> iOrder >> ittmp[0] >> ittmp[1] >> ittmp[2] >> ittmp[3] >> ittmp[4] >> ittmp[5];
    if(!mapFile) {
      Nrow[i - 1] = j;
      break;
    }
    if (iTmp > iTmp_prev){
      j = j + 1;
      iTmp_prev = iTmp;
    }else{
      Nrow[i - 1] = j;
      i = i + 1;
      j = 1;
      iTmp_prev = 1;
    }
    M[(j - 1) * 6 + i - 1] = Mtmp1;
    for(int k = 1; k <= 6; ++k){
      it[(j - 1) * 6 * 6 + (i - 1) * 6 + k - 1] = ittmp[k - 1];
    }
    //std::cout << iTmp << "\t" << Mtmp1 << "\t"  << iOrder << "\t"  << ittmp[0] << "\t"  << ittmp[1] << "\t"  << ittmp[2] << "\t"  << ittmp[3] << "\t"  << ittmp[4] << "\t"  << ittmp[5] << "\n";
  }
  mapFile.close();

  //!-----------------------------------------------------
  //! Truncate the matrix at the prescribed order Norder
  //!-----------------------------------------------------
  jmaxOrd = 0;
  for(i = 1; i <= 6; ++i){
    int kl = 0;
    for(j = 1; j <= Nrow[i - 1]; ++j){
      int sum = 0;
      for(int k = 1; k <= 6; ++k){
	sum += it[(j - 1) * 6 * 6 + (i - 1) * 6 + k - 1];
      }
      if(sum <= Norder) kl = kl+1;
      if(jmaxOrd < sum)jmaxOrd = sum;
    }
    Nrow[i - 1] = kl;
  }
}


void InputProcessing::printMap(Map *map, std::string filename){
  std::cout << "\n\nTODO:printMap!!!\n\n";
  //for(int i = 1; i <= 6 ; ++i ){
  //for(int j = 1; j <= Nrow[i - 1]; ++j){
      
  //}
  //}
}

//Reading in the order of Ncol * Npart
void InputProcessing::readICs(double *&x, std::string fileName, int &Npart, int Ncol){
  std::string filePath = paths["input"] + fileName;
  std::ifstream file(filePath.c_str());
  if (!file.is_open()){
    std::string msg = "Cannot open file " + fileName;
    Abort(msg.c_str());
  }
  double xval = 0;
  file >> Npart;
  for(int i = 0; i < Npart; ++i){
    for(int j = 0; j < Ncol; ++j){
      file >> xval;
      x[j * Npart + i] = xval;
      //x[i * Ncol + j] = xval;
    }
  }
  file.close();
}

void InputProcessing::dumpParticles(double *x, int Npart, int Ncol, int Nfreq, int iTurn, std::string ic, std::_Ios_Openmode mode, int pid){
  std::stringstream str;
  str << "PID_" << pid << ".Turn_" << iTurn;
  std::string filename = outputConfig[ic] + "." + str.str();
  std::ofstream file(filename.c_str(), mode);
  if (!file.is_open()){
    std::string msg = "Cannot open file " + outputConfig[ic];
    Abort(msg.c_str());
  }
  std::stringstream ss;
  ss.precision(16);
  ss << std::scientific;
  //ss << std::fixed;
  for(int i = 0; i < Npart; ++i){
    ss << iTurn ;
    for(int j = 0; j < Ncol; ++j){
      ss << "\t" << x[j * Npart + i];
    }
    ss << "\n";
  } 
  file << ss.str();
  file.close();
  log_in_progress = false;
}

void InputProcessing::threadFinalize(){
#if GCC_VERSION > 40800
  if (log_in_progress && log_thread.joinable()) { log_thread.join(); }
#endif
}

void InputProcessing::dumpLum(int iTurn, double Lc, double Lsl, int N,
			      double xbar_e, double ybar_e, double xbar_p, double ybar_p, double sig_x_e, double sig_y_e,
			      double sig_x_p, double sig_y_p, double *mom_x_e, double *mom_y_e, double *mom_x_p, double *mom_y_p,
			      double pxbar_e, double pybar_e, double pzbar_e, double sig_px_e, double sig_py_e, double sig_pz_e,
			      double pxbar_p, double pybar_p, double pzbar_p, double sig_px_p, double sig_py_p, double sig_pz_p, std::string sFile, std::_Ios_Openmode mode){
  std::ofstream file(outputConfig[sFile].c_str(), mode);
  if (!file.is_open()){
    std::string msg = "Cannot open file " + outputConfig[sFile];
    Abort(msg.c_str());
  }

  double sig_x = sqrt(sig_x_e * sig_x_e + sig_x_p * sig_x_p);
  double sig_y = sqrt(sig_y_e * sig_y_e + sig_y_p * sig_y_p);
  double arg = -0.50*pow(((xbar_p-xbar_e)/sig_x), 2.0) -0.50 * pow(((ybar_p-ybar_e)/sig_y), 2.0);
  double Lum = Lc * exp(arg)/(sig_x*sig_y);

  std::stringstream ss;
  ss.precision(12);
  ss << iTurn << "\t" << Lum << "\t";                //&  ! 1-2:   turn, luminosity
  ss << sig_x_e << "\t" << sig_y_e << "\t";          //&  ! 3-4:   rms size for the e-beam in x, y
  ss << sig_x_p << "\t" << sig_y_p << "\t";          //&  ! 5-6:   rms size for the p-beam in x, y
  ss << sig_x << "\t" << sig_y << "\t";              //&  ! 7-8:   combined rms size in x, y
  ss << mom_x_e[0] << "\t" << mom_x_e[1] << "\t";    //&  ! 9-10:  3rd & 4th mom. of e-beam in x
  ss << mom_y_e[0] << "\t" << mom_y_e[1] << "\t";    //&  ! 11-12: 3rd & 4th mom. of e-beam in y
  ss << mom_x_p[0] << "\t" << mom_x_p[1] << "\t";    //&  ! 13-14: 3rd & 4th mom. of p-beam in x
  ss << mom_y_p[0] << "\t" << mom_y_p[1] << "\t";    //&  ! 15-16: 3rd & 4th mom. of p-beam in y
  ss << Lsl << "\t";                                 //&  ! 17:    luminosity from slices
  ss << sig_px_e << "\t" << sig_py_e << "\t";        //&  ! 18-19: rms of px and py for the e-beam
  ss << sig_px_p << "\t" << sig_py_p << "\n";        //    ! 20-21: rms of px and py for the p-beam
  file << ss.str();
  file.close();
  
}
