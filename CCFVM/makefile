# Set the F90 Compiler
TARGETDIR="/home/kiranhegde/FlowSolvers/UNSFlow2DRun/"
#TARGETDIR="./"

FC     =gfortran
FFLAGS = -march=native  -O5 -pipe -ffree-form \
			-Wall -Wextra -Wtabs -Wintrinsics-std -Wsurprising -Waliasing \
			-Wno-unused-parameter -Wno-unused-dummy-argument \
			-Wno-unused-variable  -Warray-bounds  -fcheck=all -fbounds-check  \
                        -Waliasing -Wall -Wampersand -Warray-bounds -Wc-binding-type \
                        -Wcharacter-truncation -Wconversion -Wfunction-elimination \
                        -Wintrinsic-shadow -Wintrinsics-std \
                        -Wline-truncation -Wno-align-commons -Wno-tabs -Wreal-q-constant \
                        -Wsurprising -Wunderflow -Wunused-parameter \
                        -fbacktrace
                       # -Wimplicit-procedure -Wimplicit-interface\ 
                       # -Wuse-without-only -fbacktrace -Wrealloc-lhs -Wimplicit-interface\
                       # -fbacktrace -fdump-fortran-optimized -fdump-fortran-original
                        #-fsyntax-only -pedantic -pedantic-errors 
			#-Wrealloc-lhs-all -Wfrontend-loop-interchange -Wtarget-lifetime \


HDR = $(wildcard *.h)
SRC = $(wildcard *.F90)
OBJ = $(patsubst %.F90,%.o,$(SRC))

# Modules
MOD = data_type.mod commons.mod visc.mod inf.mod  param.mod  pri.mod   grid.mod  tools.mod
OMOD = $(patsubst %.mod,%.o,$(MOD))

TARGET  = Flo2d 

all:     $(TARGET)

data_type.mod: precision/data_type.F90
	$(FC) -c $(FFLAGS) precision/data_type.F90

%.mod: types/%.F90 data_type.mod
	$(FC) -c $(FFLAGS) types/$*.F90

%.o: %.F90 $(HDR) $(MOD) $(OMOD)
	$(FC) -c $(FFLAGS) $(INC)  $*.F90 -o $*.o

$(TARGET): $(OBJ) $(HDR) $(MOD) $(OMOD)
	$(FC) -o ${TARGETDIR}$(TARGET) $(OBJ) $(OMOD) $(LDFLAGS)

clean:
	rm -f *.o *.mod ${TARGETDIR}$(TARGET)
