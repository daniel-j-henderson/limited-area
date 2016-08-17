   program main
      use netcdf
      use mpas_file_manip
      use utils
      use params

      implicit none
      
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      ! real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      integer :: nCells, nCellsLocal, nEdges, nEdgesLocal, nVertices, nVerticesLocal, maxEdges, vertexDegree
      integer, dimension(:), pointer :: bdyMaskCell, bdyMaskEdge, bdyMaskVertex, &
                                        bdyMaskCellLocal, bdyMaskEdgeLocal, bdyMaskVertexLocal, nEdgesOnCell
      integer, dimension(:), pointer :: cell_map, icell_map, edge_map, iedge_map, vertex_map, ivertex_map
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, CellsOnVertex
      real (kind=RKIND) :: radius
      character(len=StrKIND) :: static_file = ' ', file_a = ' ', file_b = ' ', file_c = ' ', output_filename = ' ', tempstr=' ', tempstrb = ' '
      character(len=StrKIND), dimension(:), allocatable :: files, files_b, files_c
      character(len=StrKIND), dimension(NVARS_MAX_NAMELIST) :: vars_a = ' ', vars_b = ' ', vars_c = ' ', vars_r = ' '
      character(len=200) :: filename
      character(len=600) :: arg
      type(ncfile) :: ncin, ncout, nca, ncb, ncc, ncr
      integer :: i, j, l, ierr, iFile, nFiles_a=1, nFiles_b = 1, nFiles_c = 1
      integer :: aint_present, testing
      logical :: related = .false., nostatic = .false.

      namelist /filesandvariables/ static_file, file_a, vars_a, file_b, vars_b, file_c, vars_c, vars_r, output_filename     
      namelist /test/ testing 

      aint_present = access(name='namelist.limitedarea', mode=' ')      
      if (aint_present == 0) then
         open(12, file='namelist.limitedarea')
         read(12, filesandvariables)
         close(12)
      end if

      l = command_argument_count()
      i = 1
      do while (i <= l)
         call get_command_argument(i, value=tempstr)
         select case(trim(tempstr))
         case('-file_a')
            i = i+1
            call get_command_argument(i, value=file_a)
         case('-file_b')
            i = i+1
            call get_command_argument(i, value=file_b)
         case('-file_c')
            i = i+1
            call get_command_argument(i, value=file_c)
         case('-vars_a')
            i = i+1
            j = 1
            call get_command_argument(i, value=tempstrb)
            do while(tempstrb(1:1) /= '-' .and. i <= l)
               if(j <= NVARS_MAX_NAMELIST) vars_a(j) = trim(tempstrb)
               i = i+1
               j = j+1
               call get_command_argument(i, value=tempstrb)
            end do
            i = i-1 
         case('-vars_b')
            i = i+1
            j = 1
            call get_command_argument(i, value=tempstrb)
            do while(tempstrb(1:1) /= '-' .and. i <= l)
               if(j <= NVARS_MAX_NAMELIST) vars_b(j) = trim(tempstrb)
               i = i+1
               j = j+1
               call get_command_argument(i, value=tempstrb)
            end do
            i = i-1 
         case('-vars_c')
            i = i+1
            j = 1
            call get_command_argument(i, value=tempstrb)
            do while(tempstrb(1:1) /= '-' .and. i <= l)
               if(j <= NVARS_MAX_NAMELIST) vars_c(j) = trim(tempstrb)
               i = i+1
               j = j+1
               call get_command_argument(i, value=tempstrb)
            end do
            i = i-1 
         case('-vars_r')
            i = i+1
            j = 1
            call get_command_argument(i, value=tempstrb)
            do while(tempstrb(1:1) /= '-' .and. i <= l)
               if(j <= NVARS_MAX_NAMELIST) vars_r(j) = trim(tempstrb)
               i = i+1
               j = j+1
               call get_command_argument(i, value=tempstrb)
            end do
            i = i-1 
         case('-o')
            i = i+1
            call get_command_argument(i, value=output_filename)
         case('-i')
            i = i+1
            call get_command_argument(i, value=static_file)
         case('-related_files')
            related = .true.

            i = i+1
            nFiles_a = 0
            call get_command_argument(i, value=tempstrb)
            do while(tempstrb(1:1) /= '-' .and. i <= l)
               i = i+1
               nFiles_a = nFIles_a+1
               call get_command_argument(i, value=tempstrb)
            end do
            i = i - nFiles_a - 1

            allocate(files(nFiles_a))
            i = i+1
            j = 1
            call get_command_argument(i, value=tempstrb)
            do while(tempstrb(1:1) /= '-' .and. i <= l)
               files(j) = trim(tempstrb)
               i = i+1
               j = j+1
               call get_command_argument(i, value=tempstrb)
            end do
            i = i-1 
            case('-nostatic')
               nostatic = .true.
         case default
         end select
         
         i = i+1

      end do

!      write (0,*) trim(file_a)
!      write (0,*) trim(file_b)
!      write (0,*) trim(file_c)
!      do i=1, NVARS_MAX_NAMELIST
!         write (0,*) trim(vars_a(i))
!      end do
!      do i=1, NVARS_MAX_NAMELIST
!         write (0,*) trim(vars_b(i))
!      end do
!      do i=1, NVARS_MAX_NAMELIST
!         write (0,*) trim(vars_c(i))
!      end do
!      write (0,*) trim(output_filename)
!      write (0,*) trim(static_file)
      
      if (len(trim(static_file)) == 0) then
         write (0,*) "==========================================================="
         write (0,*) "   Please provide an input filename for the static file."
         write (0,*) "==========================================================="
         stop
      end if
      aint_present = access(name=static_file, mode=' ')
      if (aint_present /= 0) then
         write (0,*) "==========================================================="
         write (0,*) "      The file "//trim(static_file)//" is not present"
         write (0,*) "==========================================================="
         stop
      end if
   
      ncin%filename = trim(static_file)


      write (0,*) "Opening file "//trim(ncin%filename)//"..."
      call open_mpas_file(ncin, 'NF90_NOWRITE')

      if (len(trim(file_a)) .ne. 0) then
         nca%filename = trim(file_a)
         write (0,*) "Opening file "//trim(nca%filename)//"..."
         call open_mpas_file(nca, 'NF90_NOWRITE')
      end if

      if (len(trim(file_b)) .ne. 0) then
         ncb%filename = trim(file_b)
         write (0,*) "Opening file "//trim(ncb%filename)//"..."
         call open_mpas_file(ncb, 'NF90_NOWRITE')
      end if

      if (len(trim(file_c)) .ne. 0) then
         ncc%filename = trim(file_c)
         write (0,*) "Opening file "//trim(ncc%filename)//"..."
         call open_mpas_file(ncc, 'NF90_NOWRITE')
      end if
      
      !write (0,*) "ncids:", ncin%ncid, nca%ncid, ncout%ncid
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
      !call get_variable_1dREAL(ncin, 'xCell', xCell)
      !call get_variable_1dREAL(ncin, 'yCell', yCell)
      !call get_variable_1dREAL(ncin, 'zCell', zCell)

 
      allocate(bdyMaskCell(nCells), bdyMaskEdge(nEdges), bdyMaskVertex(nVertices), source=0)
   
      write (0,*) "Creating boundary..."
      ! call create_boundary(nCells, radius, bdyMaskCell, nEdgesOnCell, cellsOnCell, latCell, lonCell, xCell, yCell, zCell)
      call create_boundary(nCells, 1.0, bdyMaskCell, nEdgesOnCell, cellsOnCell, latCell, lonCell)

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

      do iFile=1, nFiles_a
         if(related) then
            ncout%filename = trim(output_filename)//trim(files(iFile))
            if (len(trim(output_filename)) == 0) then
               ncout%filename = 'limited-area.'//trim(files(iFile))
            else
               ncout%filename = trim(output_filename)//trim(files(iFile))
            end if
            ncr%filename = trim(files(iFile))
         else
            if (len(trim(output_filename)) == 0) then
               ncout%filename = 'limited-area.'//trim(static_file)
            else
               ncout%filename = trim(output_filename)
            end if
            ncr%filename = ' '
         end if

         write (0,*) "Creating file "//trim(ncout%filename)//"..."
         call open_mpas_file(ncout, 'CREATE')

         if (len(trim(ncr%filename)) .ne. 0) then
            write (0,*) "Opening file "//trim(ncr%filename)//"..."
            call open_mpas_file(ncr, 'NF90_NOWRITE')
         end if

         ncout%nCells = nCellsLocal
         ncout%nEdges = nEdgesLocal
         ncout%nVertices = nVerticesLocal


         call add_dimension(ncout, 'nCells', nCellsLocal)
         call add_dimension(ncout, 'nEdges', nEdgesLocal)
         call add_dimension(ncout, 'nVertices', nVerticesLocal)

         call copy_dimensions(ncin, ncout)      
         call copy_attributes(ncin, ncout)
        

         if (ncr%is_open()) then
            call copy_dimensions(ncr, ncout)
            do i=1, size(vars_r)
               if(len(trim(vars_r(i)))==0) then
                  cycle
               end if
               if (.not. ncr%contains_elem(VAR, trim(vars_r(i)))) then
                  write (0,*) "  "//trim(ncr%filename)//" is missing variable "//trim(vars_r(i))//", skipping it." 
                  cycle
               end if
               call copy_variable_defmode(ncr, ncout, vars_r(i))
            end do
         end if

         if (nca%is_open()) then
            call copy_dimensions(nca, ncout)
            do i=1, size(vars_a)
               if(len(trim(vars_a(i)))==0) then
                  cycle
               end if
               if (.not. nca%contains_elem(VAR, trim(vars_a(i)))) then
                  write (0,*) "  "//trim(nca%filename)//" is missing variable "//trim(vars_a(i))//", skipping it." 
                  cycle
               end if
               call copy_variable_defmode(nca, ncout, vars_a(i))
            end do
         end if
            
         if (ncb%is_open()) then
            call copy_dimensions(ncb, ncout)
            do i=1, size(vars_b)
               if(len(trim(vars_b(i)))==0) then
                  cycle
               end if
               if (.not. ncb%contains_elem(VAR, trim(vars_b(i)))) then
                  write (0,*) "  "//trim(ncb%filename)//" is missing variable "//trim(vars_b(i))//", skipping it." 
                  cycle
               end if
               call copy_variable_defmode(ncb, ncout, vars_b(i))
            end do
         end if
            
         if (ncc%is_open()) then
            call copy_dimensions(ncc, ncout)
            do i=1, size(vars_c)
               if(len(trim(vars_c(i)))==0) then
                  cycle
               end if
               if (.not. nca%contains_elem(VAR, trim(vars_c(i)))) then
                  write (0,*) "  "//trim(ncc%filename)//" is missing variable "//trim(vars_c(i))//", skipping it." 
                  cycle
               end if
               call copy_variable_defmode(ncc, ncout, vars_c(i))
            end do
         end if

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
      

         write (0,*) "Closing up the files..."

         ! Copy all requested variables and close the optional files
         if (ncr%is_open()) then
            do i=1, size(vars_r)
               if (.not. ncr%contains_elem(VAR, vars_r(i))) then
                  cycle
               end if
               call copy_variable_datamode(ncr, ncout, vars_r(i), cell_map, edge_map, vertex_map)
            end do
            call close_mpas_file(ncr)
         end if

         if (nca%is_open()) then
            do i=1, size(vars_a)
               if (.not. nca%contains_elem(VAR, vars_a(i))) then
                  cycle
               end if
               call copy_variable_datamode(nca, ncout, vars_a(i), cell_map, edge_map, vertex_map)
            end do
            call close_mpas_file(nca)
         end if

         if (ncb%is_open()) then
            do i=1, size(vars_b)
               if (.not. ncb%contains_elem(VAR, vars_b(i))) then
                  cycle
               end if
               call copy_variable_datamode(ncb, ncout, vars_b(i), cell_map, edge_map, vertex_map)
            end do
            call close_mpas_file(ncb)
         end if
            
         if (ncc%is_open()) then
            do i=1, size(vars_c)
               if (.not. ncc%contains_elem(VAR, vars_c(i))) then
                  cycle
               end if
               call copy_variable_datamode(ncc, ncout, vars_c(i), cell_map, edge_map, vertex_map)
            end do
            call close_mpas_file(ncc)
         end if

         call close_mpas_file(ncout)
      end do

      call close_mpas_file(ncin)
     
      write (0,*) "All Done."

   
   end program main 
