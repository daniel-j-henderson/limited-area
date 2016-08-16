
INCLUDES = $(wildcard *.o)

ifneq "$(NETCDF)" ""
        INCLUDES += -I$(NETCDF)/include
        LIBS += -L$(NETCDF)/lib
        NCLIB = -lnetcdf
        NCLIBF = -lnetcdff
        ifneq ($(wildcard $(NETCDF)/lib/libnetcdff.*), ) # CHECK FOR NETCDF4
                LIBS += $(NCLIBF)
        endif # CHECK FOR NETCDF4
        LIBS += $(NCLIB)
endif


gfortran:
		( $(MAKE) all \
        "FC = gfortran" \
        "FFLAGS = -g -fdump-core -fbacktrace -ffree-form --std=legacy -ffree-line-length-none -fdefault-real-8" \
		"LDFLAGS = ")	

ifort:
		( $(MAKE) all \
        "FC = ifort" \
        "FFLAGS = -g -debug all -traceback -autodouble" \
		"LDFLAGS = ")		
		
pgfortran:
		( $(MAKE) all \
        "FC = pgfortran" \
         "FFLAGS = " \
		"LDFLAGS = ")	



all: limited_area.o
	$(FC) $(FFLAGS) limited_area.o utils.o kd_tree.o min_heap.o file_manip.o -o limited_area $(LIBS)

limited_area.o: limited_area.f90 utils.o kd_tree.o min_heap.o file_manip.o params.o
	$(FC) $(FFLAGS) -c limited_area.f90 $(INCLUDES) 

kd_tree.o: kd_tree.f90 params.o
	$(FC) $(FFLAGS) -c kd_tree.f90 $(INCLUDES)

min_heap.o: min_heap.f90 params.o
	$(FC) $(FFLAGS) -c min_heap.f90 $(INCLUDES)

file_manip.o: file_manip.f90 utils.o params.o
	$(FC) $(FFLAGS) -c file_manip.f90 $(INCLUDES)

params.o: params.f90
	$(FC) $(FFLAGS) -c params.f90

utils.o: utils.f90 kd_tree.o min_heap.o params.o
	$(FC) $(FFLAGS) -c utils.f90 $(INCLUDES)

clean:
	rm *.o limited_area *.mod 
