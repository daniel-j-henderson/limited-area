all: limited_area.f90 utils.o file_manip.o params.o
	$(FC) $(FFLAGS) limited_area.f90 utils.o file_manip.o params.o $(INCLUDES) $(LIBS) -o ../limited_area 

file_manip.o: file_manip.f90 utils.o params.o
	$(FC) $(FFLAGS) -c file_manip.f90 $(INCLUDES)

params.o: params.f90
	$(FC) $(FFLAGS) -c params.f90

utils.o: utils.f90 params.o
	$(FC) $(FFLAGS) -c utils.f90 $(INCLUDES)

clean:
	rm *.o *.mod 
