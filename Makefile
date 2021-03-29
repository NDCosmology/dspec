#---------------------- Metal options
#OPT += -DSTELLARAGE    # SFR must also be defined. Outputs Star Age to file
#OPT += -DFEEDBACK
#OPT += -DCARBON        # Turn on the respective metal flag if gadget was run with it on.
#OPT += -DNITROGEN
#OPT += -DOXYGEN        
#OPT += -DFLORINE
#OPT += -DNEON
#OPT += -DSODIUM
#OPT += -DMAGNESIUM
#OPT += -DALUMINUM
#OPT += -DSILICON
#OPT += -DPHOSPHORUS
#OPT += -DSULFUR
#OPT += -DCHLORINE
#OPT += -DARGON
#OPT += -DPOTASSIUM
#OPT += -DCALCIUM        #Use
#OPT += -DSCANDIUM
#OPT += -DTITANIUM
#OPT += -DVANADIUM
#OPT += -DCHROMIUM       #Use
#OPT += -DMANGANESE      #Use
#OPT += -DIRON           #Use
#OPT += -DCOBALT
#OPT += -DNICKEL
#OPT += -DCOPPER
#OPT += -DZINC


#---------------------- Other options
OPT += -DPERIODIC
OPT += -DPMGRID=128


#---------------------- Cosmology options
OPT += -DLITTLE_H=0.6731
OPT += -DOMEGA_BARYON=0.0456
OPT += -DOMEGA_M0=0.315
OPT += -DOMEGA_R0=0.0
OPT += -DOMEGA_K0=0.0
OPT += -DOMEGA_DE0=0.685
OPT += -DBARYON_FRACTION=0.155

#--------------------Select Target Computer

#SYSTYPE="crc"
#SYSTYPE="ali_home"
#SYSTYPE="phillips"
SYSTYPE="jared_home"

#--------------------System Specific Settings



ifeq ($(SYSTYPE),"ali_home")
#CC       =  gcc
CC       =  mpicc
#OPTIMIZE =  -O3 -Wall -g
#OPTIMIZE =   -Wall  -gdwarf-2 -g3                    #For Debugging...
OPTIMIZE =   -Wall  -g3                               #For Debugging...
GSL_INCL =  -I/opt/local/include
GSL_LIBS =  -L/opt/local/lib  -Wl #,"-R /usr/common/pdsoft/lib"
endif




ifeq ($(SYSTYPE),"crc")
CC = mpicc
#OPTIMIZE = -Wall -g3  #Debugging
OPTIMIZE = -Wall -O3
#OPTIMIZE = -Wall -g3 -O3 -fopenmp
#OPTIMIZE = -Wall -g3 -fopenmp
GSL_INCL = -I/opt/crc/scilib/gsl/1.16/intel-14.0/include
GSL_LIBS = -L/opt/crc/scilib/gsl/1.16/intel-14.0/lib  -Wl,"-R /opt/crc/scilib/gsl/1.16/intel-14.0/lib"
#FFTW_INCL=  -I/opt/crc/scilib/fftw/2.1.5_ompi/intel/1.6.5/include
#FFTW_LIBS=  -L/opt/crc/scilib/fftw/2.1.5_ompi/intel/1.6.5/lib
#MPICHLIB =  -L/opt/crc/openmpi/1.6.5/intel-14.0/lib
MPICHLIB =  -L/opt/crc/o/openmpi/1.10.2/intel/15.0/lib

#MNI     = /opt/crc/sandbox/isuh/gadget2/minc/2.2.0
#OPT_DIR = /opt/crc/sandbox/isuh/gadget2/minc/2.2.0
#NETCDF_DIR = /opt/crc/netcdf/rhel6/4.3.0/gcc-4.4.7
#HDF_DIR    = /opt/crc/hdf/rhel6/1.8.8/gcc-4.4.6
#INCLUDE =  -I$(OPT_DIR)/include $(GSL_INCL) -I$(MNI)/include -I$(MNI)/include/volume_io     -I$(NETCDF_DIR)/include -I$(HDF_DIR)/include
#LIBS = -L$(OPT_DIR)/lib  -L$(MNI)/lib $(GSL_LIBS) $(MPICHLIB) -L$(NETCDF_DIR)/lib -L$(HD    F_DIR)/lib -lgsl -lgslcblas -lvolume_io2 -lminc2 -lnetcdf -lhdf5 #-lsrfftw_mpi -lsfftw_m    pi -lsrfftw -lsfftw
#OPT += -DNOTYPEPREFIX_FFTW

endif



ifeq ($(SYSTYPE),"jared_home")
CC = mpicc
OPTIMIZE = -Wall -g3
#OPTIMIZE = -Wall -O3
#OPTIMIZE = -Wall -O3 -fopenmp
#OPTIMIZE = -Wall -gstabs -fopenmp
GSL_INCL = -I/usr/local/include/gsl
GSL_LIBS = -L/usr/local/lib
endif



ifeq ($(SYSTYPE),"phillips")
CC = mpicc
OPTIMIZE = -Wall -O3
GSL_INCL = -I/opt/local/include
GSL_LIBS = -L/opt/local/lib -Wl
endif



#-------------------Bookkeeping
PREFIX = ./src
OBJ_DIR = ./src/obj
EXEC = dspec

OPTIONS = $(OPTIMIZE) $(OPT)
OBJS   = $(OBJ_DIR)/main.o  $(OBJ_DIR)/predict.o \
   $(OBJ_DIR)/begrun.o $(OBJ_DIR)/endrun.o  \
   $(OBJ_DIR)/init.o  $(OBJ_DIR)/io.o    \
   $(OBJ_DIR)/read_ic.o  $(OBJ_DIR)/ngb.o  \
   $(OBJ_DIR)/system.o  $(OBJ_DIR)/allocate.o   \
   $(OBJ_DIR)/forcetree.o \
   $(OBJ_DIR)/domain.o  $(OBJ_DIR)/allvars.o  \
   $(OBJ_DIR)/peano.o $(OBJ_DIR)/simplify.o \
   $(OBJ_DIR)/density.o $(OBJ_DIR)/write.o \
   $(OBJ_DIR)/error_check.o



INCL = $(PREFIX)/allvars.h $(PREFIX)/proto.h Makefile
#INCLUDE = -Iinclude -I/opt/local/include -I/usr/local/include -I/usr/include -I/usr/local/include/volume_io
INCLUDE = $(GSL_INCL)
CFLAGS = $(OPTIONS) $(GSL_INCL)

all: $(EXEC)

LIBS = $(GSL_LIBS) -lgsl -lgslcblas -lm

$(OBJ_DIR)/%.o : $(PREFIX)/%.c $(INCL)
	$(CC) $(OPTIONS) $(INCLUDE) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

#----------------Clean
clean: 
	rm -f $(OBJS) *.gch  
