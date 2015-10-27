#include "common.h"
#include "Beam.h"
#include "util/Util.h"
#include "IO.h"


class Simulator{
 public:
  InputProcessing *io;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  int num_tasks, rank, deviceId;

  Beam *beam;
  BeamParams *bParams;
  Simulator(int pid = 0);
  void Simulate(int, char **);
  bool Initialize();
  void Tracking(int &, int &);
  bool listen(int master_pid = 0);

  int* generateMapData(Map *&map, int maxLen, int Npart, int Ncol);
  MPI_Datatype getMetadataType();
  void sendMetadata(int maxLen_eqns, int maxLen, int Npart, int Ncol, BeamParams *bParams);

  void logParticles(double *h_x, int Npart, int Ncol, int iTurn);

};


struct metadata{
  int maxLen_eqns, maxLen, Npart, Ncol;
};

