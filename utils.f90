  module utils
   use params
   implicit none
   integer, parameter :: BOUNDARY1 = 2
   integer, parameter :: BOUNDARY2 = 3
   integer, parameter :: BOUNDARY3 = 4
   integer, parameter :: BOUNDARY4 = 5
   integer, parameter :: BOUNDARY5 = 6
   integer, parameter :: INSIDE = 1
   contains
    ! subroutine create_boundary(nCells, radius, bdyMaskCell, nEdgesOnCell, cellsOnCell, latCell, lonCell, xCell, yCell, zCell)
    subroutine create_boundary(nCells, radius, bdyMaskCell, nEdgesOnCell, cellsOnCell, latCell, lonCell)

      use kd_tree_mod
      use minheap_mod

      implicit none

      integer :: iCell, source_cell, v, i, j, k, npts, target_cell, inside_cell
      integer, intent(in) :: nCells
      real (kind=RKIND), dimension(:), pointer, intent(in) :: latCell, lonCell
      ! real (kind=RKIND), dimension(:), pointer, intent(in) :: xCell, yCell, zCell
      integer, dimension(:), pointer, intent(out) :: bdyMaskCell
      integer, dimension(:), pointer, intent(in) :: nEdgesOnCell
      integer, dimension(:,:), pointer, intent(in) :: cellsOnCell
      integer, dimension(:), allocatable :: boundary_cells
      real (kind=RKIND), dimension(:,:), allocatable :: cellpoints 
      real (kind=RKIND) :: mindist, dist, PI=3.141592653589793
      real (kind=RKIND), intent(in) :: radius
      real (kind=RKIND), dimension(3) :: pt

      integer, dimension(:), allocatable :: prev
      logical, dimension(:), allocatable :: unvisited
      real (kind=RKIND), dimension(:), allocatable :: distance
      type(kd_tree) :: cell_tree
      type(min_heap) :: q

      integer(kind=RKIND) :: t1, t2, t3, rate



      call system_clock(t1, rate)

      write (0,*) "   Creating cell tree..."
      ! Create cell tree
      allocate(cellpoints(3,nCells))
      do i=1,nCells
      !   cellpoints(:,i) = (/xCell(i), yCell(i), zCell(i)/)
         call con_lx(latCell(i), lonCell(i), radius, pt(1), pt(2), pt(3))
         cellpoints(:,i) = pt
      end do
      call cell_tree%create_tree(cellpoints)
      deallocate(cellpoints)

      ! Read in a list of ordered lat-lon points, do not repeat the first point at the end
      open(10, FILE='points.txt')
      read(10, *) npts
      pt = 0.0
      allocate(boundary_cells(npts))
      do i=1, npts
         read(10,*) pt(1), pt(2)
         pt = pt * PI / 180.0
         call con_lx(pt(1), pt(2), radius, pt(1), pt(2), pt(3))
         boundary_cells(i) = cell_tree%nearest_cell(pt)
      end do

      read(10,*) pt(1), pt(2)
      pt = pt * PI / 180.0
      call con_lx(pt(1), pt(2), radius, pt(1), pt(2), pt(3))
      inside_cell = cell_tree%nearest_cell(pt)

      call system_clock(t2)
      write (0,*) "   Time to create cell tree and read points:", real(t2-t1) / real(rate)

      allocate(prev(nCells), unvisited(nCells), distance(nCells))

      do i=1,npts
         ! For each path segment of adjacent boundary points...
         source_cell = boundary_cells(i)
         target_cell = boundary_cells(mod(i, npts) + 1)
         call q%create_heap(100, nCells, source_cell, real(0.0, kind=RKIND)) 
         distance = huge(dist)
         prev = 0
         unvisited = .true.
         distance(source_cell) = 0.0
         ! ...Perform Dijkstra's Algorithm to find the shortest path
         do j=1,nCells
            iCell = q%extract_min()
            if (iCell == target_cell) exit
            unvisited(iCell) = .false.
            do k = 1, nEdgesOnCell(iCell)
               v = cellsOnCell(k, iCell)
               dist = distance(iCell) + sphere_distance(latCell(iCell), lonCell(iCell), latCell(v), lonCell(v), radius)
               if (dist < distance(v)) then
                  distance(v) = dist
                  if (q%index_array(v) == 0) then
                     call q%insert(v, dist)
                  else 
                     call q%decrease_priority(v, dist)
                  end if
                  prev(v) = iCell
               end if
            end do
         end do
         iCell = target_cell
         do while(iCell .ne. source_cell)
            bdyMaskCell(iCell) = INSIDE
            iCell = prev(iCell)
         end do
         bdyMaskCell(source_cell) = INSIDE
         call q%delete_heap()
      end do
      
      call system_clock(t3)
      write (0,*) "   Time to do Dijkstra's Algorithm with the heap: ", real(t3-t2) / real(rate) 

      call system_clock(t2)
      bdyMaskCell(inside_cell) = INSIDE
      call mark_neighbors_from_source(inside_cell, INSIDE, bdyMaskCell, cellsOnCell, nEdgesOnCell)      
      call mark_neighbors_of_type(INSIDE, BOUNDARY1, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY1, BOUNDARY2, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY2, BOUNDARY3, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY3, BOUNDARY4, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY4, BOUNDARY5, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call system_clock(t3)
      write (0,*) "   Time to mark the relevant cells: ", real(t3-t2) / real(rate)


      ! Optionally make the boundary points and nearby cells a different value so they stand out in ncview, for testing purposes
      do i=1, npts
         bdyMaskCell(boundary_cells(i)) = 10
!         do j=1, nEdgesOnCell(boundary_cells(i))
!            bdyMaskCell(cellsOnCell(j, boundary_cells(i))) = 10
!         end do
      end do
      
      write (0,*) "   Time for whole find_boundary_cells routine: ", real(t3-t1) / real(rate)

   end subroutine create_boundary 

   subroutine mark_neighbors_of_type(type1, type2, mask, cellsOnCell, nEdgesOnCell)
   ! For each cell in the global mesh, if it is of type1, mark all of its
   ! unmarked neighbors as type 2
      integer, intent(in) :: type1, type2
      integer, dimension(:), pointer, intent(inout) :: mask
      integer, dimension(:), pointer, intent(in) :: nEdgesOnCell
      integer, dimension(:,:), pointer, intent(in) :: cellsOnCell

      integer :: i, v, iCell, nCells

      nCells = size(mask)

      do iCell = 1, nCells
         if (mask(iCell) .ne. type1) cycle
         do i = 1, nEdgesOnCell(iCell)
            v = cellsOnCell(i, iCell)
            if (mask(v) == 0) mask(v) = type2
         end do
      end do
   end subroutine mark_neighbors_of_type

   recursive subroutine mark_neighbors_from_source(inside_cell, tval, mask, cellsOnCell, nEdgesOnCell)
   ! Starting at some cell, expand radially, marking all unmarked cells as tval
   ! until the boundary is reached.
      integer, intent(in) :: inside_cell, tval
      integer, dimension(:), pointer, intent(inout) :: mask
      integer, dimension(:), pointer, intent(in) :: nEdgesOnCell
      integer, dimension(:,:), pointer, intent(in) :: cellsOnCell
      
      integer :: i, j, iCell
   
      do i=1, nEdgesOnCell(inside_cell)
         iCell = cellsOnCell(i, inside_cell)
         if (mask(iCell) == 0) then
            mask(iCell) = tval
            call mark_neighbors_from_source(iCell, tval, mask, cellsOnCell, nEdgesOnCell)
         end if
      end do
   end subroutine mark_neighbors_from_source

   subroutine create_local_area_map(mask, map, imap)
   ! 'map' maps the limited-area local cell (or edge, vert) to their respective global mesh
   ! source cells. 'imap' is the opposite, mapping each global mesh cell to its
   ! respective local cell, or to '0'.
      implicit none
      
      integer, dimension(:), pointer, intent(in) :: mask
      integer, dimension(:), pointer, intent(out) :: map, imap
      
      integer :: i, j

      j = 0
      do i=1, size(mask)
         if (mask(i) .ne. 0) j=j+1
      end do

      allocate(map(j), imap(size(mask)), source=0)

      j = 1
      do i=1, size(mask)
         if (mask(i) .ne. 0) then
            map(j) = i
            imap(i) = j
            j = j+1
         end if
      end do
   end subroutine create_local_area_map

   subroutine compact_field_1dINT(field, newfield, map)
   ! Use the map to create the limited-area version of a field
      integer, dimension(:), pointer :: field, newfield, map
      
      integer :: i
   
      !if (associated(newfield)) deallocate(newfield)
      !allocate(newfield(size(map)))

      do i=1, size(newfield)
         newfield(i) = field(map(i))
      end do
   end subroutine compact_field_1dINT

   subroutine compact_field_1dREAL(field, newfield, map)
      real(kind=RKIND), dimension(:), pointer :: field, newfield
      integer, dimension(:), pointer :: map
      
      integer :: i
   
      !if (associated(newfield)) deallocate(newfield)
      !allocate(newfield(size(map)))

      do i=1, size(newfield)
         newfield(i) = field(map(i))
      end do
   end subroutine compact_field_1dREAL

   subroutine compact_field_2dINT(field, newfield, map)
      integer, dimension(:), pointer, intent(in) :: map
      integer, dimension(:,:), pointer, intent(inout) :: field, newfield

      integer :: i, j
      integer, dimension(2) :: dims

      dims = shape(field)
      !if (associated(newfield)) deallocate(newfield)
      !allocate(newfield(dims(1), size(map)))

      do i=1, size(map)
      do j=1, dims(1)
         newfield(j, i) = field(j, map(i))
      end do
      end do
   end subroutine compact_field_2dINT

   subroutine compact_field_3dINT(field, newfield, map)
      integer, dimension(:), pointer :: map
      integer, dimension(:,:,:), pointer :: field, newfield

      integer :: i, j
      integer, dimension(3) :: dims

      dims = shape(field)
      !if (associated(newfield)) deallocate(newfield)
      !allocate(newfield(dims(1), dims(2), size(map)))

      do i=1, size(map)
         newfield(:,:,i) = field(:,:,map(i))
      end do
   end subroutine compact_field_3dINT

   subroutine compact_field_2dREAL(field, newfield, map)
      integer, dimension(:), pointer :: map
      real(kind=RKIND), dimension(:,:), pointer :: field, newfield

      integer :: i
      integer, dimension(2) :: dims

      dims = shape(field)
      !if (associated(newfield)) deallocate(newfield)
      !allocate(newfield(dims(1), size(map)))

      do i=1, size(map)
         newfield(:, i) = field(:, map(i))
      end do
   end subroutine compact_field_2dREAL

   subroutine compact_field_3dREAL(field, newfield, map)
      integer, dimension(:), pointer :: map
      real(kind=RKIND), dimension(:,:,:), pointer :: field, newfield

      integer :: i
      integer, dimension(3) :: dims

      dims = shape(field)
      !if (associated(newfield)) deallocate(newfield)
      !allocate(newfield(dims(1), dims(2), size(map)))

      do i=1, size(map)
         newfield(:, :, i) = field(:, :, map(i))
      end do
   end subroutine compact_field_3dREAL

   subroutine reindex_field_2dINT(field, imap)
   ! Takes an integer field whose values are element ids in the global mesh and 
   ! reindexes them to the limited-area mesh. 'field' should already be
   ! compacted.

      integer, dimension(:,:) :: field
      integer, dimension(:) :: imap

      integer :: i, j
      integer, dimension(2) :: dims

      dims = shape(field)

      do i=1, dims(2)
      do j=1, dims(1)
         if (field(j,i) > size(imap)) field(j,i) = 0 !if unused elements are indexed to nElems+1
         if (.not. field(j,i) > 0) cycle
         field(j, i) = imap(field(j, i))
      end do
      end do
   end subroutine reindex_field_2dINT


   real (kind=RKIND) function sphere_distance(lat1, lon1, lat2, lon2, radius)

      implicit none

      real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius
      real (kind=RKIND) :: arg1

      arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                 cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
      sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance

   subroutine con_lx(lat, lon, radius, x, y, z)
      implicit none

      real (kind=RKIND), intent(in) :: radius, lat, lon
      real (kind=RKIND), intent(out) :: x, y, z

      z = radius * sin(lat)
      x = radius * cos(lon) * cos(lat)
      y = radius * sin(lon) * cos(lat)
   end subroutine con_lx

end module utils
