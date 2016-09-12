   program main
      use netcdf
      use mpas_file_manip
      use utils
      use params

      implicit none
      
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:,:), allocatable :: bdy_points
      real (kind=RKIND), dimension(2) :: inside_pt
      integer :: nCells, nCellsLocal, nEdges, nEdgesLocal, nVertices, nVerticesLocal, maxEdges, vertexDegree
      integer, dimension(:), pointer :: bdyMaskCell, bdyMaskEdge, bdyMaskVertex, &
                                        bdyMaskCellLocal, bdyMaskEdgeLocal, bdyMaskVertexLocal, nEdgesOnCell
      integer, dimension(:), pointer :: cell_map, icell_map, edge_map, iedge_map, vertex_map, ivertex_map
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, CellsOnVertex
      real (kind=RKIND) :: radius
      character(len=StrKIND) :: static_file = ' '
      character(len=StrKIND) :: region_prefix
      character(len=StrKIND), dimension(:), allocatable :: input_files
      type(ncfile) :: ncin, ncout, ncr
      integer :: i, j, l, ierr, iFile, nFiles
      logical :: file_present, nostatic = .false.



      nFiles = command_argument_count()
      allocate(input_files(nFiles))
      do i=1, nFiles
         call get_command_argument(i, value=input_files(i))
         if (len(trim(static_file)) == 0) then
            if (file_contains_elem(input_files(i), VAR, 'cellsOnCell')) then
               static_file = input_files(i)
            end if
         end if
         inquire(file=trim(input_files(i)), exist=file_present)
         if (.not. file_present) then
            write (0,*) "==========================================================="
            write (0,*) "      The file "//trim(static_file)//" is not present"
            write (0,*) "==========================================================="
            stop
         end if
      end do


      if (len(trim(static_file)) == 0) then
         write (0,*) "==========================================================="
         write (0,*) "        Please provide an file which contains all" 
         write (0,*) "               the necessary static fields."
         write (0,*) "==========================================================="
         stop
      end if
      inquire(file=trim(static_file), exist=file_present)
      if (.not. file_present) then
         write (0,*) "==========================================================="
         write (0,*) "      The file "//trim(static_file)//" is not present"
         write (0,*) "==========================================================="
         stop
      end if

      
   
      ncin%filename = trim(static_file)


      write (0,*) "Opening file "//trim(ncin%filename)//"..."
      call open_mpas_file(ncin, 'NF90_NOWRITE')

      
      write (0,*) "Reading dimensions and variables..."
      call get_dimension(ncin, 'nCells', nCells)
      call get_dimension(ncin, 'nEdges', nEdges)
      call get_dimension(ncin, 'nVertices', nVertices)
      call get_dimension(ncin, 'maxEdges', maxEdges)
      call get_dimension(ncin, 'vertexDegree', vertexDegree)
      call get_attribute_REAL(ncin, 'sphere_radius', radius)

      call get_variable_1dINT(ncin, 'nEdgesOnCell', nEdgesOnCell)
      call get_variable_2dINT(ncin, 'cellsOnCell', cellsOnCell)
      call get_variable_2dINT(ncin, 'cellsOnEdge', cellsOnEdge)
      call get_variable_2dINT(ncin, 'cellsOnVertex', cellsOnVertex)
         
      call get_variable_1dREAL(ncin, 'latCell', latCell)
      call get_variable_1dREAL(ncin, 'lonCell', lonCell)
 
      allocate(bdyMaskCell(nCells), bdyMaskEdge(nEdges), bdyMaskVertex(nVertices), source=0)
   
      write (0,*) "Creating boundary..."
      call open_pointfile('points.txt', bdy_points, inside_pt, region_prefix)
      call create_boundary(nCells, 1.0, bdy_points, inside_pt, bdyMaskCell, nEdgesOnCell, cellsOnCell, latCell, lonCell)

      do i=1, nEdges
         bdyMaskEdge(i) = max(bdyMaskCell(cellsOnEdge(1, i)), bdyMaskCell(cellsOnEdge(2, i)))
      end do

      do i=1, nVertices
         do j=1, vertexDegree
            if (bdyMaskCell(cellsOnVertex(j,i)) > bdyMaskVertex(i)) bdyMaskVertex(i) = bdyMaskCell(cellsOnVertex(j,i))
         end do
      end do

      call create_local_area_map(bdyMaskCell, cell_map, icell_map)
      call create_local_area_map(bdyMaskEdge, edge_map, iedge_map)
      call create_local_area_map(bdyMaskVertex, vertex_map, ivertex_map)
   
      nCellsLocal = size(cell_map)
      nEdgesLocal = size(edge_map)
      nVerticesLocal = size(vertex_map)

      write (0,*) "nCells, edges, vertices", ncin%nCells, ncin%nEdges, ncin%nVertices
      write (0,*) "nCellsLocal, nEdgesLocal, nVerticesLocal:", nCellsLocal, nEdgesLocal, nVerticesLocal

      do iFile=1, nFiles
         ncout%filename = trim(region_prefix)//'.'//trim(input_files(iFile))
         ncr%filename = trim(input_files(iFile))

         write (0,*) "Creating file "//trim(ncout%filename)//"..."
         call open_mpas_file(ncout, 'CREATE')
         if (len(trim(ncr%filename)) .ne. 0) then
            if (trim(ncr%filename) .ne. trim(ncin%filename)) then
               write (0,*) "Opening file "//trim(ncr%filename)//"..."
               call open_mpas_file(ncr, 'NF90_NOWRITE')
            else
               ncr = ncin
            end if
            if (ncr%nCells .ne. -1 .and. ncr%nCells .ne. ncin%nCells) then ! the file has nCells but not the same nCells as the static
               write (0,*) "---------------------------------------------------------------------------------"
               write (0,*) '  The file '//trim(ncr%filename)//' you provided has different dimensions '
               write (0,*) '  appears to be on a different grid than the static file '//trim(static_file)
               write (0,*) '  Please make sure all files are on the same grid'
               write (0,*) "---------------------------------------------------------------------------------"
               stop
            end if
            if (ncr%nEdges .ne. -1 .and. ncr%nEdges .ne. ncin%nEdges) then
               write (0,*) "---------------------------------------------------------------------------------"
               write (0,*) '  The file '//trim(ncr%filename)//' you provided has different dimensions '
               write (0,*) '  appears to be on a different grid than the static file '//trim(static_file)
               write (0,*) '  Please make sure all files are on the same grid'
               write (0,*) "---------------------------------------------------------------------------------"
               stop
            end if
            if (ncr%nVertices .ne. -1 .and. ncr%nVertices .ne. ncin%nVertices) then
               write (0,*) "---------------------------------------------------------------------------------"
               write (0,*) '  The file '//trim(ncr%filename)//' you provided has different dimensions '
               write (0,*) '  appears to be on a different grid than the static file '//trim(static_file)
               write (0,*) '  Please make sure all files are on the same grid'
               write (0,*) "---------------------------------------------------------------------------------"
               stop
            end if
         end if

         ncout%nCells = nCellsLocal
         ncout%nEdges = nEdgesLocal
         ncout%nVertices = nVerticesLocal


         call add_dimension(ncout, 'nCells', nCellsLocal)
         call add_dimension(ncout, 'nEdges', nEdgesLocal)
         call add_dimension(ncout, 'nVertices', nVerticesLocal)

         !call copy_dimensions(ncr, ncout)      
         call copy_attributes(ncr, ncout)
        
         ! Copy all non-static variables into the new file
         if (ncr%is_open()) then
            call copy_dimensions(ncr, ncout)
            do i=1, ncr%nvars
               if(len(trim(ncr%vars(i)))==0 .or. ncr%vars(i)(1:7) == 'indexTo' .or. is_static(trim(ncr%vars(i)))) then
                  cycle
               end if
               call copy_variable_defmode(ncr, ncout, ncr%vars(i))
            end do
         end if

         ! Definition of a file having 'static' fields is that it contains
         ! cellsOnCell
         if (ncr%contains_elem(VAR, 'cellsOnCell')) then
            nostatic = .false.
         else
            nostatic = .true.
         end if

         ! If a file has static fields, then copy all static fields into the new
         ! limited-area file also. Otherwise, just end define mode and move on.
         if (.not. nostatic) then

            allocate(bdyMaskCellLocal(nCellsLocal)) 
            call compact_field_1dINT(bdyMaskCell, bdyMaskCellLocal, cell_map)
            call create_variable_1dINT(ncout, 'bdyMaskCell', 'nCells')

            allocate(bdyMaskEdgeLocal(nEdgesLocal))
            call compact_field_1dINT(bdyMaskEdge, bdyMaskEdgeLocal, edge_map)
            call create_variable_1dINT(ncout, 'bdyMaskEdge', 'nEdges')

            allocate(bdyMaskVertexLocal(nVerticesLocal))
            call compact_field_1dINT(bdyMaskVertex, bdyMaskVertexLocal, vertex_map)
            call create_variable_1dINT(ncout, 'bdyMaskVertex', 'nVertices')
            ! these 3 variables are just simpler to make separately
            call create_variable_1dINT(ncout, 'indexToCellID', 'nCells')
            call create_variable_1dINT(ncout, 'indexToEdgeID', 'nEdges')
            call create_variable_1dINT(ncout, 'indexToVertexID', 'nVertices')

            call copyandcompact_static_fields(ncin, ncout, cell_map, edge_map, vertex_map, icell_map, iedge_map, ivertex_map)

            call put_variable_1dINT(ncout, bdyMaskCellLocal, 'bdyMaskCell')
            call put_variable_1dINT(ncout, bdyMaskEdgeLocal, 'bdyMaskEdge')
            call put_variable_1dINT(ncout, bdyMaskVertexLocal, 'bdyMaskVertex')

            bdyMaskCellLocal = (/(i, i=1,ncout%nCells)/)
            bdyMaskEdgeLocal = (/(i, i=1,ncout%nEdges)/)
            bdyMaskVertexLocal = (/(i, i=1,ncout%nVertices)/)

            call put_variable_1dINT(ncout, bdyMaskCellLocal, 'indexToCellID')
            call put_variable_1dINT(ncout, bdyMaskEdgeLocal, 'indexToEdgeID')
            call put_variable_1dINT(ncout, bdyMaskVertexLocal, 'indexToVertexID')
         else
            ierr = nf90_enddef(ncout%ncid)
            if (ierr /= NF90_NOERR) then
               write(0,*) '*********************************************************************************'
               write(0,*) 'Error ending define mode'
               write(0,*) 'ierr = ', ierr
               write(0,*) '*********************************************************************************'
            end if
         end if
      

         write (0,*) "Closing up the new file..."

         ! Copy all variables from the input file to the output
         if (ncr%is_open()) then
            do i=1, ncr%nvars
               if (.not. is_static(trim(ncr%vars(i)))) &
                  call copy_variable_datamode(ncr, ncout, ncr%vars(i), cell_map, edge_map, vertex_map)
            end do
            call close_mpas_file(ncr)
         end if


         call close_mpas_file(ncout)
      end do

      !call close_mpas_file(ncin)
     
      write (0,*) "All Done."

   
   end program main 
