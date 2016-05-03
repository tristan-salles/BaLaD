# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2014
#
# Description : Makefile for Steady State
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TARGET = InitPhase
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES = Interpolation.f90 Morphometrics.f90 TIN_computation.f90 \
	Triangular_network.f90 Facc_algorithm.f90 SteadyState.f90
		
FOBJS=$(FSOURCES:.f90=.o)
 
OBJS=./*.o 
  
.PHONY: dist module dust clobber

dist : module 

module : 	
	@echo
	@echo "*************************************************"
	@echo	
	@echo "Build modules and library for Steady State."
	@echo	
	@echo "*************************************************"
	make $(TARGETLIB)
	
$(TARGETLIB) : $(FOBJS)
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done
				
%.o : %.f90
	$(F90) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) ${KDTREEFLAGS} -c $< -o $@ 
	
dust :
	rm -f *.bak *~ *.o *.mod *.original 

clobber : dust