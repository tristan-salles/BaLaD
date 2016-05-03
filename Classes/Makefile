# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2014 
#
# Description : Makefile for Classes
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

FTRIA = ./Triangle

TARGET = Classes
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES= MPI_parameters.f90 File_parameters.f90 Mesh_parameters.f90 \
	Flow_parameters.f90 Time_parameters.f90 Forces_parameters.f90 

FOBJS=$(FSOURCES:.f90=.o)

CSOURCES= $(FTRIA)/TriangleCall.c

COBJS=$(CSOURCES:.c=.o)
OBJS=./*.o $(FTRIA)/*.o 

.PHONY: dist module dust clobber 

dist : module

module :	
	@echo "*************************************************"
	@echo	
	@echo "Build modules and library for Classes."
	@echo	
	@echo "*************************************************"
	make $(TARGETLIB)
	
$(TARGETLIB): $(COBJS) $(FOBJS)
	cp *.o $(OBJDIR)
	cp $(FTRIA)/*.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done

%.o : %.f90
	$(F90) $(FFLAGS) ${KDTREEFLAGS} -c $< -o $@ 
	
%.o : %.c
	$(CC) $(CFLAGS) -I$(TRIANGLELIB) -c $< -o $@   
	
dust:
	rm -fv *.bak *~ *.o *.mod $(FTRIA)/*.o *.original 
	
clobber: dust
