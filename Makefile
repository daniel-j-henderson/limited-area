
INCLUDES = $(wildcard *.o)

LIBS = $(shell nc-config --libs)

INCLUDES += -I$(shell nc-config --includedir)

all:
	@echo "*********************************************"
	@echo "  Use a target, such as ifort or gfortran"
	@echo "*********************************************"

gfortran:
	cd src; $(MAKE) FC="gfortran" \
	INCLUDES="$(INCLUDES)" \
	LIBS="$(LIBS)" \
	FFLAGS="-ffree-form --std=legacy -ffree-line-length-none -fdefault-real-8" 

ifort:
	cd src; $(MAKE) FC="ifort" \
	INCLUDES="$(INCLUDES)" \
	LIBS="$(LIBS)" \
	FFLAGS="-g -traceback -autodouble" 
						
clean:
	cd src; $(MAKE) clean
	rm limited_area
