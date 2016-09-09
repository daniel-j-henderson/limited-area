# LIMITED-AREA
**A tool to generate limited area MPAS meshes.**

##BUILD
Just make the makefile with your preferred compiler as the target (I support ifort or gfortran right now), e.g. 'make gfortran'.

##RUN
1. Run the utility on the command line. Provide as command line arguments all of the netcdf files you'd like to have a regional copy of. At least one of these files must contain the static fields for that mesh.

2. Provide a list of points that will make up the boundary of your limited area mesh
  - the first line of the file should be the name of the region (it will be the filename prefix for all the output files)
  - these points should be in degrees
  - the points must be ordered, as the utility will simply 'connect the dots' using the points provided
  - the second line of the file should be the number of boundary points provided
  - after the list of boundary points, you must provide a point which is guaranteed to be inside the boundary (if it is outside, you'll get the complement of your desired region)
  a sample is shown in the Examples directory
