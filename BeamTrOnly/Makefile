EXECUTABLE      := exec
NODES		:= 2
_OBJS		:= Beam.o InputProcessing.o util/Util.o 
CFILES		:= Simulator.cu
OPTIONS		:= -std=c++11 #-g -G #--ptxas-options=-v #-g -G
ODIR 		:= obj
SDIR 		:= .
CMINPACK_SDIR	:= cminpack
CMINPACK_ODIR	:= obj/cminpack
MPI_SDIR	:= mpi
MPI_ODIR	:= obj/mpi
INC 		:= 
CC 		:= nvcc
OBJS		:= $(patsubst %,$(ODIR)/%,$(_OBJS))
CCFLAGS 	:= -arch=sm_35 -lmpi  #-lgomp -Xcompiler -fopenmp 
INCLUDE := -I /usr/local/include 
LIBS    := -L /usr/local/lib


$(ODIR)/%.o:$(SDIR)/%.cpp 
	 $(CC)  ${OPTIONS} -c -o $@ $< $(CFLAGS)	

all : ${OBJS} 
	$(CC) -o ${EXECUTABLE} ${INCLUDE} ${LIBS} ${OPTIONS} ${OBJS} $(CCFLAGS) ${CFILES} 	
clean:
	rm ${OBJS}
	rm output/*
run: 
	mpirun -np ${NODES} ${EXECUTABLE}
