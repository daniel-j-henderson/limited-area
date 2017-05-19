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
      integer :: i, j, l, ierr, iFile, nFiles, region_type
      logical :: file_present, nostatic = .false.
      real(kind=RKIND), dimension(:), allocatable :: rparams



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
      call close_mpas_file(ncin)

      allocate(bdyMaskCell(nCells), bdyMaskEdge(nEdges), bdyMaskVertex(nVertices), source=UNMARKED)
   
      write (0,*) "Creating boundary..."
      call open_pointfile('points.txt', bdy_points, inside_pt, region_prefix, region_type, rparams)
      write (0,*) trim(region_prefix), region_type, rparams
      
      if(region_type == RCUSTOM) then
         call create_boundary_custom(nCells, 1.0_RKIND, bdy_points, inside_pt, bdyMaskCell, nEdgesOnCell, cellsOnCell, latCell, lonCell)
      else
         call create_boundary_from_region(bdyMaskCell, region_type, rparams, &
                                       latCell, lonCell, radius, cellsOnCell, nEdgesOnCell) 
      end if

      deallocate(nEdgesOnCell)
      deallocate(cellsOnCell)

      do i=1, nEdges
         if (bdyMaskCell(cellsOnEdge(1,i)) == UNMARKED) then
            bdyMaskEdge(i) = bdyMaskCell(cellsOnEdge(2,i))
         else if (bdyMaskCell(cellsOnEdge(2,i)) == UNMARKED) then
            bdyMaskEdge(i) = bdyMaskCell(cellsOnEdge(1,i))
         else
            bdyMaskEdge(i) = min(bdyMaskCell(cellsOnEdge(2,i)), bdyMaskCell(cellsOnEdge(1,i)))
         end if
      end do

      deallocate(cellsOnEdge)

      do i=1, nVertices
         do j=1, vertexDegree
            if (bdyMaskCell(cellsOnVertex(j,i)) == UNMARKED) then
               cycle
            else if (bdyMaskVertex(i) == UNMARKED) then
               bdyMaskVertex(i) = bdyMaskCell(cellsOnVertex(j,i)) 
            end if
            if (bdyMaskCell(cellsOnVertex(j,i)) < bdyMaskVertex(i)) bdyMaskVertex(i) = bdyMaskCell(cellsOnVertex(j,i))
         end do
      end do

      deallocate(cellsOnVertex)

      call create_local_area_map(bdyMaskCell, cell_map, icell_map)
      call create_local_area_map(bdyMaskEdge, edge_map, iedge_map)
      call create_local_area_map(bdyMaskVertex, vertex_map, ivertex_map)
   
      nCellsLocal = size(cell_map)
      nEdgesLocal = size(edge_map)
      nVerticesLocal = size(vertex_map)

      do iFile=1, nFiles
         ncout%filename = trim(region_prefix)//'.'//trim(input_files(iFile))
         ncin%filename = trim(input_files(iFile))

         write (0,*) "Creating file "//trim(ncout%filename)//"..."
         call open_mpas_file(ncout, 'CREATE')
         if (len(trim(ncin%filename)) .ne. 0) then
            write (0,*) "Opening file "//trim(ncin%filename)//"..."
            call open_mpas_file(ncin, 'NF90_NOWRITE')
         end if
         if (ncin%nCells .ne. -1 .and. ncin%nCells .ne. size(icell_map)) then ! the file has nCells but not the same nCells as the static
            write (0,*) "---------------------------------------------------------------------------------"
            write (0,*) '  The file '//trim(ncin%filename)//' you provided has different dimensions '
            write (0,*) '  appears to be on a different grid than the static file '//trim(static_file)
            write (0,*) '  Please make sure all files are on the same grid'
            write (0,*) "---------------------------------------------------------------------------------"
            stop
         end if
         if (ncin%nEdges .ne. -1 .and. ncin%nEdges .ne. size(iedge_map)) then
            write (0,*) "---------------------------------------------------------------------------------"
            write (0,*) '  The file '//trim(ncin%filename)//' you provided has different dimensions '
            write (0,*) '  appears to be on a different grid than the static file '//trim(static_file)
            write (0,*) '  Please make sure all files are on the same grid'
            write (0,*) "---------------------------------------------------------------------------------"
            stop
         end if
         if (ncin%nVertices .ne. -1 .and. ncin%nVertices .ne. size(ivertex_map)) then
            write (0,*) "---------------------------------------------------------------------------------"
            write (0,*) '  The file '//trim(ncin%filename)//' you provided has different dimensions '
            write (0,*) '  appears to be on a different grid than the static file '//trim(static_file)
            write (0,*) '  Please make sure all files are on the same grid'
            write (0,*) "---------------------------------------------------------------------------------"
            stop
         end if

         ncout%nCells = nCellsLocal
         ncout%nEdges = nEdgesLocal
         ncout%nVertices = nVerticesLocal


         call add_dimension(ncout, 'nCells', nCellsLocal)
         call add_dimension(ncout, 'nEdges', nEdgesLocal)
         call add_dimension(ncout, 'nVertices', nVerticesLocal)

         !call copy_dimensions(ncr, ncout)      
         call copy_attributes(ncin, ncout)
        
         ! Copy all non-static variables into the new file
         call copy_dimensions(ncin, ncout)
         do i=1, ncin%nvars
            if(len(trim(ncin%vars(i)))==0 .or. ncin%vars(i)(1:7) == 'indexTo' .or. is_static(trim(ncin%vars(i)))) then
               cycle
            end if
            call copy_variable_defmode(ncin, ncout, ncin%vars(i))
         end do

         ! Definition of a file having 'static' fields is that it contains
         ! cellsOnCell
         if (ncin%contains_elem(VAR, 'cellsOnCell')) then
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
            if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_enddef', .true., 'main', ncout%filename)
         end if
      
         deallocate(icell_map)
         deallocate(iedge_map)
         deallocate(ivertex_map)

         write (0,*) "Closing up the new file..."

         ! Copy all variables from the input file to the output
         do i=1, ncin%nvars
            if (trim(ncin%vars(i)) == 'indexToCellID') cycle
            if (trim(ncin%vars(i)) == 'indexToVertexID') cycle
            if (trim(ncin%vars(i)) == 'indexToEdgeID') cycle
            if (.not. is_static(trim(ncin%vars(i)))) &
               call copy_variable_datamode(ncin, ncout, ncin%vars(i), cell_map, edge_map, vertex_map)
         end do
         call close_mpas_file(ncin)


         call close_mpas_file(ncout)
      end do

      write (0,*) "All Done."

   
   end program main 
