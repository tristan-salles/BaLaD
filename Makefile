# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2014 
# Description : Makefile for BADLANDS
# ============================================================================
TOP=$(shell pwd)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

DIRMODS= Classes Visualisation Initialisation SteadyState \
	LayerStatus Flows ProcessDriver

SOURCES = badlands.f90
OBJS=$(SOURCES:.f90=.o)

.PHONY : all dist plugin dust clobber

all: dist

dist: 
	@echo
	@echo "*************************************************"
	@echo "BADLANDS Author Tristan Salles "
	@echo "*************************************************"
	@echo
	@mkdir -p $(BUILDDIR)	
	@mkdir -p $(OBJDIR)
	@mkdir -p $(MODDIR)
	@mkdir -p $(LIBDIR)
	@mkdir -p bin
	for i in $(DIRMODS) ; do   \
    	  ( cd $$i ; make dist) ;       \
	done
	@echo "*************************************************"
	@echo	
	@echo "Build BADLANDS binary."
	@echo	
	@echo "*************************************************"
	@$(if $(wildcard badlands.o),rm -f badlands.o,)	
	make $(EXEC)

plugin : 
	cd $(DIRPLUG); make plugin; 
	@echo "*************************************************"
	@echo	
	@echo "BADLANDS shared library created."
	@echo	
	@echo "*************************************************"

$(EXEC) :	$(OBJS)
	$(F90) $(FFLAGS)  $(FOXFLAGS) $(H5FLAGS) ${KDTREEFLAGS} -I$(TRIANGLELIB) -o $@ $^ \
	$(TRIANGLELIB)/*.o $(LDFLAGS) -lBADLANDS  $(H5LDFLAGS) $(H5LIBS) $(LDFOXFLAGS) ${KDTREELDFLAGS} $(KDTREELIBS) 
	@echo "*************************************************"
	@echo	
	@echo "BADLANDS updated in ./bin/."
	@echo	
	@echo "*************************************************"

%.o : %.f90
	$(F90) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) ${KDTREEFLAGS} -c $< -o $@ 
	$(AR) $(LIBDIR)/libBADLANDS.a $(OBJDIR)/*.o
	
dust :
	for i in $(DIRMODS) ; do   \
    	( cd $$i ; make dust) ;       \
	done
	$(foreach module,$(MODULES),cd $(module); cd - ; )
	rm -fv *~ *.bak *.o *.mod *.original

clobber : dust
	for i in $(DIRMODS) ; do   \
    	( cd $$i ; make clobber) ;   \
	done	
	rm -rfv $(BUILDDIR)
