# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2014
#
# Description : Makefile for Initialisation
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TARGET = InitPhase
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES = Read_Strata.f90 Read_Hydro.f90 Read_Time.f90 \
	Read_Controls.f90 Init_Porosity.f90 Init_Hemipelagic.f90 \
	IsostaticFlexure.f90 CheckPoint.f90 Init_model.f90
	
FOBJS=$(FSOURCES:.f90=.o)
 
OBJS=./*.o 
  
.PHONY: dist module dust clobber

dist : module 

module : 	
	@echo
	@echo "*************************************************"
	@echo	
	@echo "Build modules and library for Initialisation"
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