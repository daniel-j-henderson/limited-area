# LIMITED-AREA
**A tool to generate limited area MPAS meshes.**

##BUILD
Just make the makefile with your preferred compiler as the target (I support ifort or gfortran right now), e.g. 'make gfortran'.

##RUN
1. Run the utility on the command line. Provide as command line arguments all of the netcdf files you'd like to have a regional copy of. At least one of these files must contain the static fields for that mesh.

2. Provide a points.txt file in the run directory
  
  a) Use a functionally defined region 
    - the first line of the file should be the name of your region (it will be the filename prefix for all the output files)
    - the second line should be the type of region you'd like. Presently, 'circle' and 'ellipse' are supported
    - the third line should be a list of the necessary parameters for the chosen region
      'circle' parameters: center_lat, center_lon, radius
      'ellipse' parameters: center_lat, center_lon, a_lat, a_lon, semiminor_axis_angle
        note: 'a' is the point on one end of the major axis. The semiminor axis angle is essentially the length of the semiminor axis divided by the radius of the MPAS mesh. All angles should be in degrees.
  
  b) Create a custom region from your own pointset
    - the first line of the file should be the name of the region (it will be the filename prefix for all the output files)
    - the second line of the file should be the number of boundary points provided
    - then, list the boundary points
    - these points should be in degrees
    - the points must be ordered, as the utility will simply 'connect the dots' using the points provided
    - after the list of boundary points, you must provide a point which is guaranteed to be inside the boundary (if it is outside, you'll get the complement of your desired region)
    
  A sample of each type of pointlist is shown in the Examples directory
