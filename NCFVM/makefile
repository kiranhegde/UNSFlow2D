SHELL=/bin/sh
#
# This makefile was produced by /usr/bin/ftnmgen
# at 05:07:18 PM on 02/29/00
#
# If it is invoked by the command line
#	make -f makefile
# it will compile the fortran modules indicated by SRCS into the object
# modules indicated by OBJS and produce an executable named mis.
#
# To remove all the objects but leave the executables use the command line
#	make -f makefile clean
#
# To remove everything but the source files use the command line
#	make -f makefile clobber
#
# To remove the source files created by /usr/bin/ftnmgen
# and this makefile use the command line
#	make -f makefile void
#
# The parameters SRCS and OBJS should not need to be changed.  If, however,
# you need to add a new source file add the name of the source file to the
# SRCS parameter and add the name of the resulting object file to the OBJS
# parameter.
#

#SRCS= All_modules.f90 average.f90 geometric.f90  common.f90 gradient.f90 timestep.f90 farfield_flux.f90 fvresidual.f90   main.f90  lusgs.f90 read_input.f90 solid_flux.f90 vanleer_flux.f90 roe_flux.f90 limiter.f90 ausm_up_flux.f90	
#OBJS= All_modules.o  average.o geometric.o common.o gradient.o timestep.o farfield_flux.o fvresidual.o  main.o  lusgs.o read_input.o  solid_flux.o vanleer_flux.o roe_flux.o limiter.o ausm_up_flux.o

SRCS= All_modules.f90 geometric.f90 gradient.f90 limiter.f90 common.f90 timestep.f90 farfield_flux.f90 fvresidual.f90  lusgs.f90  main.f90  read_input.f90 solid_flux.f90 vanleer_flux.f90 
OBJS= All_modules.o   geometric.o gradient.o limiter.o common.o timestep.o farfield_flux.o fvresidual.o  lusgs.o main.o  read_input.o  solid_flux.o vanleer_flux.o

# Tunable parameters
#
# FC		Name of the fortran compiling system to use
# FFLAGS	Flags to the compiler
#
#  Origin 200: -64 -extend_source -mips4 -r10000 -O3
#  O2        : -extend_source -mips4 -r5000 -O3
#  Challenge : -64 -extend_source -mips4 -SWP:=ON -O3 -OPT:fold_arith_limit=2599
#
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
#

#FC  =		f95
<<<<<<< HEAD
FC  =	gfortran	
=======
FC  =		ifort
>>>>>>> a688afb70276d1c80f99d1548c4045d94d229af4
FFLAGS =   -O5
#FFLAGS =  -C
#FFLAGS =  -openmp  -parallel
#FFLAGS =  -ipo -O3
#FFLAGS =   -g -C
#FFLAGS =  -e -C
#FFLAGS = -e -fast -pad -xipo -stackvar
#FFLAGS =  -ipo 
#FFLAGS = -extend_source -O3 
#LDFLAGS =	
#LIBS =		
CMD =		Flo2d

# Rules to build mis.
#
all:		$(CMD)
$(CMD):		$(OBJS)
	   $(FC) $(LDFLAGS) $(FFLAGS) -o $(CMD) $(OBJS) $(LIBS)

All_modules.o:     All_modules.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS)  All_modules.f90 $(LIBS)

main.o:     main.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) main.f90 $(LIBS)

#ausm_up_flux.o: ausm_up_flux.f90
#	   $(FC) -c $(LDFLAGS) $(FFLAGS) ausm_up_flux.f90 $(LIBS)

#average.o:     average.f90
#	   $(FC) -c $(LDFLAGS) $(FFLAGS) average.f90 $(LIBS)

common.o:   common.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) common.f90 $(LIBS)

farfield_flux.o:  farfield_flux.f90 
	   $(FC) -c $(LDFLAGS) $(FFLAGS) farfield_flux.f90 $(LIBS)

fvresidual.o:   fvresidual.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) fvresidual.f90 $(LIBS)

geometric.o:   geometric.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) geometric.f90 $(LIBS)

gradient.o:   gradient.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) gradient.f90 $(LIBS)

timestep.o:   timestep.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) timestep.f90 $(LIBS)

#hcusp_flux.o:   hcusp_flux.f90
#	   $(FC) -c $(LDFLAGS) $(FFLAGS) hcusp_flux.f90 $(LIBS)

#kfvs_flux.o:   kfvs_flux.f90
#	   $(FC) -c $(LDFLAGS) $(FFLAGS) kfvs_flux.f90 $(LIBS)

limiter.o:   limiter.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) limiter.f90 $(LIBS)

lusgs.o:  lusgs.f90 
	   $(FC) -c $(LDFLAGS) $(FFLAGS) lusgs.f90 $(LIBS)

read_input.o:  read_input.f90 
	   $(FC) -c $(LDFLAGS) $(FFLAGS) read_input.f90 $(LIBS)

#result.o:  result.f90 
#	   $(FC) -c $(LDFLAGS) $(FFLAGS) result.f90 $(LIBS)

#roe_flux.o:  roe_flux.f90 
#	   $(FC) -c $(LDFLAGS) $(FFLAGS) roe_flux.f90 $(LIBS)

#smooth.o:  smooth.f90 
#	   $(FC) -c $(LDFLAGS) $(FFLAGS) smooth.f90 $(LIBS)

solid_flux.o:  solid_flux.f90 
	   $(FC) -c $(LDFLAGS) $(FFLAGS) solid_flux.f90 $(LIBS)

vanleer_flux.o:  vanleer_flux.f90 
	   $(FC) -c $(LDFLAGS) $(FFLAGS) vanleer_flux.f90 $(LIBS)

#visc.o:  visc.f90 
#	   $(FC) -c $(LDFLAGS) $(FFLAGS) visc.f90 $(LIBS)

clean:
	-rm -f $(OBJS) *.mod

clobber:	clean
	-rm -f $(CMD)

void:	clobber
	-rm -f $(SRCS) makefile
