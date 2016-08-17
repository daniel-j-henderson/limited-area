# LIMITED-AREA
**A tool to generate limited area MPAS meshes.**

##BUILD
Just make the makefile with your preferred compiler as the target (I support ifort or gfortran right now), e.g. 'make gfortran'.

##RUN
1. Edit the namelist called namelist.limitedarea
  - specify the .nc file in which all the static variables can be found
  - specify the output filename
  - optionally, provide up to 20 variables from up to 3 different files that you'd like to have in your output
  and/or
  - use command line arguments, eg. ./limited_area -i static.nc -file_a file.nc -vars_a ter pressure -o output.nc
  - use command line argument -nostatic to skip the output of static fields

  RELATED MODE: If you'd like to take a lot of related files on the same grid (for instance, a few dozen update files or boundary condition files) and make a limited-area version of each without re-running the program from scratch each time, simply use the -related_files <*> argument on the command line. This may not be specified in the namelist. This will create a limited-area version of each of the files in the related_files argument, each including the variables specified in the vars_r namelist/command line option. The output files will be named with the -o or output_filename option field appended to the input file name.  

2. Provide a list of points that will make up the boundary of your limited area mesh
  - these points should be in degrees
  - the points must be ordered, as the utility will simply 'connect the dots' using the points provided
  - the first line of the file should be the number of boundary points provided
  - after the list of boundary points, you must provide a point which is guaranteed to be inside the boundary (if it is outside, you'll get the complement of your desired region)
  a sample is shown in the Examples directory
