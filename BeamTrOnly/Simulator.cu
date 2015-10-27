#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "Simulator.h"

#include <cuda_runtime.h>
#include "util/cudaArchUtil.h"
#include "util/cudaDebugUtil.h"
#include "util/cudaUtil.h"
#include "util/cudaTimerUtil.h"
#include "util/cudaMemoryUtil.h"

#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>


#include "GPUBeam/Map.cu"

#define TAG  0

#define BUFSIZE 256
enum MPI_XCHANGE_TYPE{TRACKING_ONLY, EXIT, PROBE};
static std::map<MPI_XCHANGE_TYPE, std::string> MPI_DATA;

Simulator::Simulator(int pid){
  this->rank = 0;
  this->deviceId = 0;
  this->num_tasks = 1;
  io = new InputProcessing;
  bParams = new BeamParams;
  beam = new Beam;
  MPI_DATA[EXIT]="EXIT";
  MPI_DATA[PROBE]="PROBE";
  MPI_DATA[TRACKING_ONLY]="TRACKING_ONLY";  
}

bool Simulator::Initialize(){
  beam->initialize(bParams);
  return true;
}

//!--------------------------------------------------------------
//! Generate the map from the simulation parameters if Mef == 1.
//!--------------------------------------------------------------
void Simulator::Simulate(int argc, char **argv){
  quad::timer::event_pair timer_node;


  if(!Initialize()){
    Abort("Initialize Error");
  }

  //! create the linear map for e-beam if Mef == 1
  if(bParams->Mef == 1){
    bParams->Norder_e = 1;           //! linear map must be of order 1
    Abort("genLinMap not implemented!!!\n");
  }


  quad::timer::start_timer(&timer_node);
  
  //GF equations from file
  io->ReadMap(bParams->gfEqns_file, beam->fmap->Nrow, beam->fmap->M, beam->fmap->it, bParams->Norder_e, bParams->jmaxord_e);
  int maxLen_eqns = Util::maxval<int>(beam->fmap->Nrow, 6);
  
  beam->allocateMap(beam->eqns, maxLen_eqns);
  memcpy(beam->eqns->Nrow, beam->fmap->Nrow, 6 * sizeof(int));
  
  for(int j = 1; j <= maxLen_eqns; ++j){
    for(int i = 1; i <= 6; ++i){
      beam->eqns->M[(j - 1) * 6 + i - 1] = beam->fmap->M[(j - 1) * 6 + i - 1];
      for(int k = 1; k <= 6; ++k){
	beam->eqns->it[(j - 1) * 6 * 6 + (i - 1) * 6 + k -1] = beam->fmap->it[(j - 1) * 6 * 6 + (i - 1) * 6 + k - 1];
      }
    }
  }
  
  io->ReadMap(bParams->Me_file, beam->fmap->Nrow, beam->fmap->M, beam->fmap->it, bParams->Norder_e, bParams->jmaxord_e);  
  double *x_p = beam->x_p;
  int maxLen_p = 0;
  int maxLen_e = Util::maxval<int>(beam->fmap->Nrow, 6);
  beam->allocateMap(beam->map_e, maxLen_e);
  memcpy(beam->map_e->Nrow, beam->fmap->Nrow, 6 * sizeof(int));
  for(int j = 1; j <= maxLen_e; ++j){
    for(int i = 1; i <= 6; ++i){
      beam->map_e->M[(j - 1) * 6 + i - 1] = beam->fmap->M[(j - 1) * 6 + i - 1];
      for(int k = 1; k <= 6; ++k){
	beam->map_e->it[(j - 1) * 6 * 6 + (i - 1) * 6 + k -1] = beam->fmap->it[(j - 1) * 6 * 6 + (i - 1) * 6 + k - 1];
      }
    }
  }

  //!---------------------------------------------
  //! Print the map for the e-beam to the screen.
  //!---------------------------------------------
  io->printMap(beam->map_e, "out");
  quad::timer::stop_timer(&timer_node, "ReadMap+MapInit");

  //!----------------------------------------------------------
  //! Read in the ICs for the e-beam from ICe_file (ICef <> 1)
  //! or generate it from the parameters read in (ICef == 1).
  //!----------------------------------------------------------
  double *x_e = beam->x_e;
  
  quad::timer::start_timer(&timer_node);
  if(bParams->ICef != 1){
    io->readICs(x_e, bParams->ICe_file, bParams->Npart_e, NCOL);
  }else{
    beam->genICs(bParams, ELECTRON);
  }

  io->dumpParticles(x_e, bParams->Npart_e, NCOL, bParams->Nfreq, 0, "IC_e");

  quad::timer::stop_timer(&timer_node, "ICs");


  if(bParams->iTrackOnly != 1){
    //! create the linear map for p-beam if Mpf == 1
    if(bParams->Mpf == 1){
      bParams->Norder_p = 1;           //! linear map must be of order 1
      Abort("genLinMap not implemented!!!\n");
    } 

    //!-------------------------------------------------------------------
    //! Read in the map from matrix_file and truncate it to the minimum 
    //! of Norder and the jmaxord (maximum order of the map) if Mpf <> 1.
    //!-------------------------------------------------------------------
    io->ReadMap(bParams->Mp_file, beam->fmap->Nrow, beam->fmap->M, beam->fmap->it, bParams->Norder_p, bParams->jmaxord_p);
    maxLen_p = Util::maxval<int>(beam->fmap->Nrow, 6);
    beam->allocateMap(beam->map_p, maxLen_p);
    memcpy(beam->map_p->Nrow, beam->fmap->Nrow, 6 * sizeof(int));

    for(int j = 1; j <= maxLen_p; ++j){
      for(int i = 1; i <= 6; ++i){
	beam->map_p->M[(j - 1) * 6 + i - 1] = beam->fmap->M[(j - 1) * 6 + i - 1];
	for(int k = 1; k <= 6; ++k){
	  beam->map_p->it[(j - 1) * 6 * 6 + (i - 1) * 6 + k -1] = beam->fmap->it[(j - 1) * 6 * 6 + (i - 1) * 6 + k - 1];
	}
      }
    }

    //!------------------------------------------------------------
    //! Read in the ICs for the p-beam from ICp_file (ICpf <> 1)
    //! or generate it from the parameters read in (ICpf == 1).
    //!------------------------------------------------------------
    
    if(bParams->ICpf != 1){
      io->readICs(x_p, bParams->ICp_file, bParams->Npart_p, NCOL);
    }else{
      beam->genICs(bParams, PROTON);
    }
    
    io->dumpParticles(x_p, bParams->Npart_p, NCOL, bParams->Nfreq, 0, "IC_p");
  }

  std::stringstream ss;
  if(bParams->iTrackOnly != 1){
    ss << "================================================================\n";
    ss << "                   BEAM-BEAM SIMULATION \n";
    ss << "================================================================\n";    
  }else{
    ss << "================================================================\n";
    ss << "                   TRACKING SIMULATION \n";
    ss << "================================================================\n";    
  }

  ss << "------------------------------------------------------\n";
  ss << "                  Simulation parameters\n";
  ss << "------------------------------------------------------\n";
  ss << "Number of map iterations       : " << bParams->Niter << "\n";
  ss << "Freqency of printing out       : " << bParams->Nfreq << "\n";

  if(bParams->iTrackOnly != 1){
    ss << "Number of slices per beam      : " << bParams->N << "\n";
    ss << "Beam-beam effect model         : " << bParams->iRegime << "\n";
    ss << "------------------------------------------------------\n";
    if (bParams->iRegime == 1) {
      ss << "                       WEAK-STRONG\n";
    }else if (bParams->iRegime == 2) {
      ss << "                       STRONG-STRONG\n";
    }else if (bParams->iRegime == 3) {
      ss << "                       NO COLLISION: compute rms\n";
    }else{
      ss << "                       NO COLLISION\n";
    }
    ss << "------------------------------------------------------\n";
  }

  ss << "\n------------------------------------------------------\n";
  ss << "                     Beam parameters\n";
  ss << "------------------------------------------------------\n";

  ss << "---------\n" ;
  ss << " e-beam \n";
  ss << "--------- \n" ;
  ss << "   Energy                      : " << bParams->E_e << "eV\n";
  ss << "   Matrix file                 : " << bParams->Me_file << "\n";
  ss << "   IC file                     : " << bParams->ICe_file << "\n";
  ss << "   Number of test particles    : " << bParams->Npart_e << "\n";
  ss << "   Map up to order             : " << bParams->jmaxord_e << "\n";
  ss << "   Requested map order         : " << bParams->Norder_e << "\n";
  ss << "   Horizontal rms size         : " << bParams->sig_x0_e << "m\n";
  ss << "   Vertical rms size           : " << bParams->sig_y0_e << "m\n";
  ss << "   Number of particles in beam : " << bParams->N_e << "\n";
  ss << "================================================================\n";

  if(bParams->iTrackOnly != 1){
    ss << "--------- \n";
    ss << " p-beam \n";
    ss << "--------- \n";
    ss << "   Energy                      : " << bParams->E_p << "eV\n";
    ss << "   Matrix file                 : " << bParams->Mp_file << "\n";
    ss << "   IC file                     : " << bParams->ICp_file << "\n";
    ss << "   Number of test particles    : " << bParams->Npart_p << "\n";
    ss << "   Map up to order             : " << bParams->jmaxord_p << "\n";
    ss << "   Requested map order         : " << bParams->Norder_p << "\n";
    ss << "   Horizontal rms size         : " << bParams->sig_x0_p << "m\n";
    ss << "   Vertical rms size           : " << bParams->sig_y0_p << "m\n";
    ss << "   Number of particles in beam : " << bParams->N_p << "\n";
  }
  ss << "================================================================\n";
  //std::cout << ss.str() << "\n"; 


  //!************************
  //! TRACKING ONLY (1 beam) 
  //!************************
  if(bParams->iTrackOnly == 1){
    char buff[BUFSIZE];
    for(int node_id = 0; node_id < num_tasks; ++node_id){
      if(node_id != 0){
	strncpy(buff, MPI_DATA[TRACKING_ONLY].c_str(), BUFSIZE);
	MPI_Send(buff, BUFSIZE, MPI_CHAR, node_id, TAG, MPI_COMM_WORLD);
      }
    }      

    if(bParams->isGPU){
      std::cout << "Executing on GPU!\n";
      double *d_x;
      int *d_Nrow = 0, *d_itM = 0;
      int *d_eqnsNrow = 0, *d_eqnsitM = 0;

      int *h_eqnsitM = generateMapData(beam->eqns, maxLen_eqns, bParams->Npart_e, NCOL);
      int *h_itM = generateMapData(beam->map_e, maxLen_e, bParams->Npart_e, NCOL);

      sendMetadata(maxLen_eqns, maxLen_e, bParams->Npart_e, NCOL, bParams);
      //Send Eqns and Map data
      for(int node_id = 0; node_id < num_tasks; ++node_id){
	if(node_id != 0){
	  MPI_Send(h_eqnsitM, maxLen_eqns * NCOL * (NCOL + 2), MPI_INT, node_id, TAG, MPI_COMM_WORLD);
	  MPI_Send(h_itM, maxLen_e * NCOL * (NCOL + 2), MPI_INT, node_id, TAG, MPI_COMM_WORLD);
	  MPI_Send(beam->eqns->Nrow, NCOL, MPI_INT, node_id, TAG, MPI_COMM_WORLD);
	  MPI_Send(beam->map_e->Nrow, NCOL, MPI_INT, node_id, TAG, MPI_COMM_WORLD);
	}
      }

      //Partition particles between nodes
      double *h_x = 0;
      int Npart = 0;
      int numParticlesPerNode = bParams->Npart_e / num_tasks;
      for(int node_id = 0; node_id < num_tasks; ++node_id){
	int offset = node_id * numParticlesPerNode;
	int count = numParticlesPerNode;
	if(node_id == num_tasks - 1){
	  count = bParams->Npart_e - offset;
	}

	double *hx = new double[count * NCOL];

	for(int ii = 0; ii < NCOL; ++ii){
	  memcpy(hx + ii * count, x_e + ii * bParams->Npart_e + offset, sizeof(double) * count);
	}
	if(node_id != 0){
	  MPI_Send(hx, count * NCOL, MPI_DOUBLE, node_id, TAG, MPI_COMM_WORLD);
	}else{
	  h_x = hx;
	  Npart = count;
	}
      }
      initDeviceMemory(d_eqnsitM, h_eqnsitM,
		       d_eqnsNrow, beam->eqns->Nrow,
		       d_x, h_x,
		       maxLen_eqns, Npart, NCOL);

      initDeviceMemory(d_itM, h_itM,
		       d_Nrow, beam->map_e->Nrow,
		       d_x, h_x,
		       maxLen_e, Npart, NCOL);

      double *dOpx = 0;
      QuadDebug(cudaMalloc((void **)&dOpx, sizeof(double) * Npart * NCOL));
      int *dOutOfBound = 0;
      QuadDebug(cudaMalloc((void **)&dOutOfBound, sizeof(int) * Npart));
      thrust::device_ptr<int> dev_ptr(dOutOfBound);    
      thrust::fill(dev_ptr, dev_ptr + Npart, (int) 0);

      quad::timer::event_pair timer0;
      quad::timer::start_timer(&timer0);
    
      double time = 0;
      for(int iTurn = 1; iTurn <= bParams->Niter; ++iTurn){	
	double exec_time = applyMapGPU(dOutOfBound,
				       d_itM, d_eqnsitM, 
				       d_x, dOpx,
				       beam->map_e->Nrow, beam->eqns->Nrow, 
				       h_x, 
				       maxLen_e, maxLen_eqns, 
				       Npart, NCOL, bParams, iTurn);
	time += exec_time;
	dumpBeamByThread(this, d_x, h_x, Npart, NCOL, iTurn, "dump.ebeam", std::ios::app);
      }
      io->threadFinalize();
      quad::timer::stop_timer(&timer0, "GPU Tracking");

      ss.str("");
      ss << "Tracking took " << time/bParams->Niter << " ms  per turn in " << hostname << " (Rank = " << rank << ")\n" ;
      std::cout << ss.str();

      cudaFree(dOpx);

      //io->dumpParticles(h_x, Npart, NCOL, bParams->Nfreq, bParams->Niter, "dump.ebeam");
    }else{
      std::cout << "Executing on CPU!\n";
      quad::timer::event_pair timer0;
      quad::timer::start_timer(&timer0);
    
      int iTurn = 1;
      double *xi = new double[NCOL * bParams->Npart_e];
      for(iTurn = 1; iTurn <= bParams->Niter; ++iTurn){
	memcpy(xi, x_e, NCOL * bParams->Npart_e * sizeof(double));
	beam->applyM(beam->map_e->M, beam->map_e->it, x_e, beam->map_e->Nrow, maxLen_e, bParams->Npart_e, bParams->Norder_e);
      
	//TODO - xi to be updated every iteration
	if(bParams->isSymTr && bParams->NSympFreq > 0 && (iTurn % bParams->NSympFreq == 0)){
	  beam->newtonIter(beam->eqns->M, beam->eqns->it, xi, x_e, beam->eqns->Nrow, maxLen_eqns, bParams->Npart_e, bParams->Norder_e);
	}

	if(bParams->Nfreq > 0 && (iTurn % bParams->Nfreq) == 0){
	  std::stringstream ss;
	  ss << std::setprecision(16);
	  ss << iTurn << " turns finished\n";
	  std::cout << ss.str() << "\n";
	  io->dumpParticles(x_e, bParams->Npart_e, NCOL, bParams->Nfreq, iTurn, "dump.ebeam");	
	}
      }
      quad::timer::stop_timer(&timer0, "CPU applyMap");
      io->dumpParticles(x_e, bParams->Npart_e, NCOL, bParams->Nfreq, bParams->Niter, "dump.ebeam");	
    }
  }else{
    //!*************************
    //! TRACKING AND COLLISION
    //!*************************    

    quad::timer::start_timer(&timer_node);
    Tracking(maxLen_e, maxLen_p);
    quad::timer::stop_timer(&timer_node, "Tracking");
  }
}



void Simulator::Tracking(int &maxLen_e, int &maxLen_p){
}  


int* 
Simulator::generateMapData(Map *&map, int maxLen, int Npart, int Ncol){
  double *M = new double[maxLen * Ncol];
  int *it = new int[maxLen * Ncol * Ncol];
  for(int i = 0; i < maxLen; ++i){
    for(int j = 0; j < Ncol; ++j){
      M[j * maxLen + i] = map->M[i * 6 + j];
      for(int k = 0; k < 6; ++k)
	it[j * maxLen * 6 + i * 6 + k] = map->it[i * 6 * 6 + j * 6 + k];

    }
  }

  int *_it_ = new int[maxLen * Ncol * (Ncol + 2)];//additional 2 (8-B) for M  
  for(int j = 1; j <= maxLen; ++j){
    for(int i = 1; i <= Ncol; ++i){
      for(int k = 1; k <= Ncol; ++k){
	//it[(i - 1) * maxLen * Ncol + (j - 1) * Ncol + k - 1] = map->it[(j - 1) * 6 * 6 + (i - 1) * 6 + k -1];
	_it_[(i - 1) * maxLen * (Ncol + 2) + (j - 1) * (Ncol + 2) + k - 1] = map->it[(j - 1) * 6 * 6 + (i - 1) * 6 + k -1];
	//std::cout << map->it[(j - 1) * 6 * 6 + (i - 1) * 6 + k -1] << "\n";
      }
      int *tmp = (int *)&map->M[(j - 1) * Ncol + i - 1];
      _it_[(i - 1) * maxLen * (Ncol + 2) + (j - 1) * (Ncol + 2) + 6] = tmp[0];
      _it_[(i - 1) * maxLen * (Ncol + 2) + (j - 1) * (Ncol + 2) + 7] = tmp[1];
      //std::cout << map->M[(j - 1) * Ncol + i - 1] << "\n";
    }
  }

  delete it,M;
  return _it_;
}


bool Simulator::listen(int master_pid){
  bool listenStatus = false;
  char buff[BUFSIZE];
  MPI_Status stat;
  MPI_Recv(buff, BUFSIZE, MPI_CHAR, master_pid, TAG, MPI_COMM_WORLD, &stat);
  std::stringstream ss;
  double start = MPI_Wtime();
  if(strncmp(buff, MPI_DATA[EXIT].c_str(), MPI_DATA[EXIT].length()) == 0){
    //Exit
    listenStatus = false;
  }else if(strncmp(buff, MPI_DATA[TRACKING_ONLY].c_str(), MPI_DATA[TRACKING_ONLY].length()) == 0){
    MPI_Status status;
    metadata meta_data;
    MPI_Datatype metadata_type = getMetadataType();
    MPI_Recv(&meta_data, 1, metadata_type, master_pid, META_TAG, MPI_COMM_WORLD, &status);


    
    int *h_eqnsitM = new int[meta_data.maxLen_eqns * meta_data.Ncol * (meta_data.Ncol + 2 )];
    int *h_itM = new int[meta_data.maxLen * meta_data.Ncol * (meta_data.Ncol + 2 )];
    int *h_eqnsNrow = new int[meta_data.Ncol];
    int *h_Nrow = new int[meta_data.Ncol];
    
    MPI_Recv(h_eqnsitM, meta_data.maxLen_eqns * meta_data.Ncol * (meta_data.Ncol + 2 ), MPI_INT, master_pid, TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(h_itM, meta_data.maxLen * meta_data.Ncol * (meta_data.Ncol + 2 ), MPI_INT, master_pid, TAG, MPI_COMM_WORLD, &status);

    MPI_Recv(h_eqnsNrow, meta_data.Ncol, MPI_INT, master_pid, TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(h_Nrow, meta_data.Ncol, MPI_INT, master_pid, TAG, MPI_COMM_WORLD, &status);



    //Receive particle info
    // Probe for an incoming message from process zero
    MPI_Probe(master_pid, TAG, MPI_COMM_WORLD, &status);

    int size = 0;
    MPI_Get_count(&status, MPI_DOUBLE, &size);

    double *h_x = new double[size];
    MPI_Recv(h_x, size, MPI_DOUBLE, master_pid, TAG, MPI_COMM_WORLD, &status);

    double *d_x;
    int *d_Nrow = 0, *d_itM = 0;
    int *d_eqnsNrow = 0, *d_eqnsitM = 0;

    int Npart = size/meta_data.Ncol;
    initDeviceMemory(d_eqnsitM, h_eqnsitM,
		     d_eqnsNrow, h_eqnsNrow,
		     d_x, h_x,
		     meta_data.maxLen_eqns, Npart, meta_data.Ncol);
    
    initDeviceMemory(d_itM, h_itM,
		     d_Nrow, h_Nrow,
		     d_x, h_x,
		     meta_data.maxLen, Npart, meta_data.Ncol);

    double *dOpx = 0;
    QuadDebug(cudaMalloc((void **)&dOpx, sizeof(double) * Npart * meta_data.Ncol));
    int *dOutOfBound = 0;
    QuadDebug(cudaMalloc((void **)&dOutOfBound, sizeof(int) * Npart));
    thrust::device_ptr<int> dev_ptr(dOutOfBound);    
    thrust::fill(dev_ptr, dev_ptr + Npart, (int) 0);
  
    quad::timer::event_pair timer_node;
    quad::timer::start_timer(&timer_node);
    
    double time = 0;
    for(int iTurn = 1; iTurn <= bParams->Niter; ++iTurn){      
      double exec_time = applyMapGPU(dOutOfBound,
				     d_itM, d_eqnsitM, 
				     d_x, dOpx,
				     h_Nrow, h_eqnsNrow, 
				     h_x, 
				     meta_data.maxLen, meta_data.maxLen_eqns, 
				     Npart, meta_data.Ncol, bParams, iTurn);
      time += exec_time;
      
      dumpBeamByThread(this, d_x, h_x, Npart, meta_data.Ncol, iTurn, "dump.ebeam", std::ios::app);
    }
    io->threadFinalize();
          
    quad::timer::stop_timer(&timer_node, "GPU Tracking");

    ss.str("");
    ss << "Tracking took " << time/bParams->Niter << " ms per turn  in " << hostname << " (Rank = " << rank << ")\n" ;
    std::cout << ss.str();

    //write to rank file
    //io->dumpParticles(h_x, Npart, NCOL, bParams->Nfreq, bParams->Niter, "dump.ebeam");
    cudaFree(dOpx);

    listenStatus = true;
  }
  return listenStatus;
}

MPI_Datatype 
Simulator::getMetadataType(){
  MPI_Datatype metadata_type;
  
  const int num_ele = 1;
  int blocklens[num_ele];
  MPI_Aint disp[num_ele];
  MPI_Datatype types[num_ele];

  blocklens[0] = 4;
  disp[0] = 0;
  types[0] = MPI_INT;

  MPI_Type_create_struct( num_ele, blocklens, disp, types, &metadata_type );
  //MPI_Type_struct( num_ele, blocklens, disp, types, &metadata_type );
  MPI_Type_commit(&metadata_type);
  return metadata_type;
}


void 
Simulator::sendMetadata(int maxLen_eqns, int maxLen, int Npart, int Ncol, BeamParams *bParams){
  metadata data;
  data.maxLen_eqns = maxLen_eqns;data.maxLen = maxLen;data.Npart = Npart;data.Ncol = Ncol;
  
  MPI_Datatype metadata_type = getMetadataType();
  

  for(int node_id = 0; node_id < num_tasks; ++node_id){
    if(node_id != 0){
      MPI_Send(&data, 1, metadata_type, node_id, META_TAG, MPI_COMM_WORLD);
    }
  }
}


int main(int argc, char **argv){
  int len = 0, rc = 0;
  std::stringstream ss;

  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    std::stringstream ss;
    ss << "Error starting MPI program. Terminating.\n";
    PrintlnPID(ss.str(), 0);
    MPI_Abort(MPI_COMM_WORLD, rc);
  }


  Simulator sim;
  MPI_Comm_size(MPI_COMM_WORLD,&sim.num_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&sim.rank);
  MPI_Get_processor_name(sim.hostname, &len);

  std::map<std::string, int> nodeMap;
  for(int node_id = 0; node_id < sim.num_tasks; ++node_id){
    char* buf = new char[MPI_MAX_PROCESSOR_NAME];
    memcpy(buf, sim.hostname, sizeof(char) * MPI_MAX_PROCESSOR_NAME);
    MPI_Bcast(buf, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, node_id, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(nodeMap.find(buf) == nodeMap.end())
      nodeMap[buf] = 0;
    else
      nodeMap[buf]++;
    if(node_id == sim.rank){
      sim.deviceId = nodeMap[buf];
      QuadDebug(cudaSetDevice(sim.deviceId));
      //ss << node_id << " " << buf << " " << nodeMap[buf] << "\n";
    }
  }
  quad::timer::event_pair timer_node;
  quad::timer::start_timer(&timer_node);
  sim.io->ReadInput(sim.bParams);
  quad::timer::stop_timer(&timer_node, "ReadInput");

  //CPU Execution & MPI == STOP
  if(!sim.bParams->isGPU && sim.num_tasks > 1){
    std::stringstream ss;
    ss << "Error: [CPU Execution does not run on MPI Cluster]\n";
    ss << "Usage: ./exec\n";
    std::cout << ss.str();
    MPI_Finalize();
    return 0;
  }
  
  if(sim.bParams->isGPU && sim.bParams->iTrackOnly == 1 && (sim.bParams->Npart_e % (ILP * sim.num_tasks)) != 0){
    std::stringstream ss;
    ss << "Error: GPU implementation requires #. of particles in the simulation should be a multiple of ILP * #. of GPUs  ";
    ss << ILP << " * " << sim.num_tasks << "\n";
    std::cout << ss.str();
    MPI_Finalize();
    return 0;
  }

  if(sim.bParams->iTrackOnly == 1 && ILP !=4  && sim.bParams->isSymTr){
    std::stringstream ss;
    ss << "Error: Symplectic Tracking is not implemented with ILP = " << ILP << "\n";
    std::cout << ss.str();
    MPI_Finalize();
    return 0;
  }
  

  if(sim.rank == 0){
    IO::FileCreatePID(0);
    sim.Simulate(argc, argv);
    char buff[BUFSIZE];
    for(int node_id = 0; node_id < sim.num_tasks; ++node_id){
      if(node_id != 0){
	strncpy(buff, MPI_DATA[EXIT].c_str(), BUFSIZE);
	MPI_Send(buff, BUFSIZE, MPI_CHAR, node_id, TAG, MPI_COMM_WORLD);
      }
    }  
  }else{

    while(true){
      bool active = sim.listen();
      if(!active)break;
    }
  }
  MPI_Finalize();
  return 0;
}
