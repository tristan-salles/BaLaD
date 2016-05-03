# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2014
#
# Description : Makefile for Flows Computation
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TARGET = FlowWalkers
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES =  ODEsolver.f90 FWsFunctions.f90 FWsEroDep.f90 OverlandFWs.f90 \
	 RegolithFormation.f90 SedTransfert.f90 HillSlope.f90 MassWasting.f90 LandscapeEvol.f90 \
	 ChannelFWs.f90

OBJS=$(FSOURCES:.f90=.o)
  
.PHONY: dist module dust clobber

dist : module 

module : 	
	@echo
	@echo "*************************************************"
	@echo	
	@echo "Build modules and library for Flow Walkers"
	@echo	
	@echo "*************************************************"
	make $(TARGETLIB)
	
$(TARGETLIB) : $(OBJS)
	#$(AR) $(TARGETLIB) $(OBJS) 
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done
				
%.o : %.f90
	$(F90) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) ${KDTREEFLAGS} -c $< -o $@ 
	
dust :
	rm -f *.o *.mod *.original

clobber : dust