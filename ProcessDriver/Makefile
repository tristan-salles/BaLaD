# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2014
#
# Description : Makefile for Process Driver
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TARGET = Process
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES = UDWplugin.f90 RunProcesses.f90
	
FOBJS=$(FSOURCES:.f90=.o)
 
OBJS=./*.o 
  
.PHONY: dist module dust clobber

dist : module 

module : 	
	@echo
	@echo "*************************************************"
	@echo	
	@echo "Build modules and library for Process Driver."
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
	rm -f *.bak  *.o *.mod *.original 

clobber : dust