  module utils
   use params
   implicit none
   integer, parameter :: BOUNDARY1 = 1
   integer, parameter :: BOUNDARY2 = 2
   integer, parameter :: BOUNDARY3 = 3
   integer, parameter :: BOUNDARY4 = 4
   integer, parameter :: BOUNDARY5 = 5
   integer, parameter :: BOUNDARY6 = 6
   integer, parameter :: BOUNDARY7 = 7
   integer, parameter :: INSIDE = 0
   integer, parameter :: UNMARKED = -1
   
   contains
    subroutine create_boundary(nCells, radius, bdy_pts, inside_pt, bdyMaskCell, nEdgesOnCell, cellsOnCell, latCell, lonCell)

      use minheap_mod

      implicit none

      integer :: iCell, source_cell, v, i, j, k, npts, target_cell, inside_cell
      integer, intent(in) :: nCells
      real (kind=RKIND), dimension(:), pointer, intent(in) :: latCell, lonCell
      real (kind=RKIND), dimension(:,:), intent(in) :: bdy_pts
      real (kind=RKIND), dimension(2) :: inside_pt
      ! real (kind=RKIND), dimension(:), pointer, intent(in) :: xCell, yCell, zCell
      integer, dimension(:), pointer, intent(out) :: bdyMaskCell
      integer, dimension(:), pointer, intent(in) :: nEdgesOnCell
      integer, dimension(:,:), pointer, intent(in) :: cellsOnCell
      integer, dimension(:), allocatable :: boundary_cells
      real (kind=RKIND), dimension(:,:), allocatable :: cellpoints 
      real (kind=RKIND) :: mindist, dist, angle, minangle 
      real (kind=RKIND), intent(in) :: radius
      real (kind=RKIND), dimension(3) :: pt, pta, ptb, v1, v2, v3

      integer, dimension(:), allocatable :: prev
      logical, dimension(:), allocatable :: unvisited
      real (kind=RKIND), dimension(:), allocatable :: distance
      type(min_heap) :: q

      integer(kind=RKIND) :: t1, t2, t3, rate
      real(kind=RKIND) :: temp
      

      call system_clock(t1, rate)

      ! Interpolate the boundary points provided onto the MPAS mesh
      npts = size(bdy_pts(1,:))
      allocate(boundary_cells(npts))
      do i=1, size(bdy_pts(1,:))
         boundary_cells(i) = nearest_cell_path(bdy_pts(1,i), bdy_pts(2,i), &
                  merge(boundary_cells(i-1), 1, i > 1), nCells, 10, nEdgesOnCell, cellsOnCell, latCell, lonCell)
      end do
      inside_cell = nearest_cell_path(inside_pt(1), inside_pt(2), &
                  boundary_cells(1), nCells, 10, nEdgesOnCell, cellsOnCell, latCell, lonCell)
         
      call system_clock(t2)
      bdyMaskCell = UNMARKED
      allocate(prev(nCells), unvisited(nCells), distance(nCells))

      ! Follow-the-line Algorithm :: A greedy algorithm that is greedy on angle
      ! between a cell and the great-circle arc from source to target 
      do i=1, npts
         source_cell = boundary_cells(i) ! beginning of segment of boundary
         target_cell = boundary_cells(mod(i, npts) + 1) ! end of segment of boundary
         call con_lx(latCell(source_cell), lonCell(source_cell), 1.0, pta(1), pta(2), pta(3))
         call con_lx(latCell(target_cell), lonCell(target_cell), 1.0, ptb(1), ptb(2), ptb(3))
         pta = cross(pta, ptb)
         temp = mag(pta)
         pta = pta / temp ! now pta = unit normal vector to the plane containing the arc from source to target
         iCell = source_cell
         do while(iCell /= target_cell) 
            bdyMaskCell(iCell) = INSIDE 
            minangle = huge(1.0) !These angles will be angles between the
                                 !cell center and the great-circle arc from source to target
            mindist = sphere_distance(latCell(iCell), lonCell(iCell),&
                                      latCell(target_cell), lonCell(target_cell), radius)
            do j=1, nEdgesOnCell(iCell)
               v = cellsOnCell(j, iCell) ! v = the jth neighbor of iCell
               dist = sphere_distance(latCell(v), lonCell(v), &
                                      latCell(target_cell), lonCell(target_cell), radius)
               if (dist > mindist) cycle ! if v was further away than iCell, skip it, regardless of its angle
                                         ! with the arc from source to target
               call con_lx(latCell(v), lonCell(v), 1.0, pt(1), pt(2), pt(3))
               angle = dot(pta, pt) ! both pt and pta are radius 1.0 so (pta dot pt) = cos(theta), where theta is angle between normal vector and pt
               angle = abs(PI / 2 - acos(angle))
               if (angle < minangle) then ! find iCell's neighbor that is both nearer the target cell than iCell 
                                          ! and minimizes the angle between the other neighbors and the arc from source to target
                  minangle = angle
                  k = v
               end if
            end do
            iCell = k
         end do
      end do


      ! Greedy Algorithm :: Alternative to the Dijkstra Algorithm
!      do i=1, npts
!         source_cell = boundary_cells(i)
!         target_cell = boundary_cells(mod(i, npts) + 1)
!         iCell = source_cell
!         do while(iCell /= target_cell) 
!            bdyMaskCell(iCell) = INSIDE
!            mindist = huge(1.0)
!            angle = sphere_distance(latCell(iCell), lonCell(iCell), latCell(target_cell), lonCell(target_cell), radius)
!            do j=1, nEdgesOnCell(iCell)
!               v = cellsOnCell(j, iCell)
!               dist = sphere_distance(latCell(v), lonCell(v), &
!                                      latCell(target_cell), lonCell(target_cell), radius)
!               if (dist > angle) cycle
!               if (dist < mindist) then
!                  mindist = dist
!                  k = v
!               end if
!            end do
!            iCell = k
!         end do
!      end do


      ! Dijkstra's Algorithm :: Produces correct boundaries, but they can
      ! sometimes cut into the desired area and don't necessarily follow the
      ! 'straight line' from source to destination
!      do i=1,npts
!         ! For each path segment of adjacent boundary points...
!         source_cell = boundary_cells(i)
!         target_cell = boundary_cells(mod(i, npts) + 1)
!         call q%create_heap(100, nCells, source_cell, real(0.0, kind=RKIND)) 
!         distance = huge(dist)
!         prev = 0
!         unvisited = .true.
!         distance(source_cell) = 0.0
!         ! ...Perform Dijkstra's Algorithm to find the shortest path
!         do j=1,nCells
!            iCell = q%extract_min()
!            if (iCell == target_cell) exit
!            unvisited(iCell) = .false.
!            do k = 1, nEdgesOnCell(iCell)
!               v = cellsOnCell(k, iCell)
!               dist = distance(iCell) + sphere_distance(latCell(iCell), lonCell(iCell), latCell(v), lonCell(v), radius)
!               if (dist < distance(v)) then
!                  distance(v) = dist
!                  if (q%index_array(v) == 0) then
!                     call q%insert(v, dist)
!                  else 
!                     call q%decrease_priority(v, dist)
!                  end if
!                  prev(v) = iCell
!               end if
!            end do
!         end do
!         iCell = target_cell
!         do while(iCell .ne. source_cell)
!            bdyMaskCell(iCell) = INSIDE
!            iCell = prev(iCell)
!         end do
!         bdyMaskCell(source_cell) = INSIDE
!         call q%delete_heap()
!      end do
      

      bdyMaskCell(inside_cell) = INSIDE 
      call mark_neighbors_from_source(inside_cell, INSIDE, bdyMaskCell, cellsOnCell, nEdgesOnCell)      
      call mark_neighbors_of_type(INSIDE, BOUNDARY1, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY1, BOUNDARY2, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY2, BOUNDARY3, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY3, BOUNDARY4, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY4, BOUNDARY5, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY5, BOUNDARY6, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY6, BOUNDARY7, bdyMaskCell, cellsOnCell, nEdgesOnCell)

      call system_clock(t3)

      ! Optionally make the boundary points and nearby cells a different value so they stand out in ncview, for testing purposes
!      do i=1, npts
!         bdyMaskCell(boundary_cells(i)) = 10
!         do j=1, nEdgesOnCell(boundary_cells(i))
!            bdyMaskCell(cellsOnCell(j, boundary_cells(i))) = 10
!         end do
!      end do
      
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
            if (mask(v) == UNMARKED) mask(v) = type2
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
         if (mask(iCell) == UNMARKED) then
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
         if (mask(i) .ne. UNMARKED) j=j+1
      end do

      allocate(map(j), imap(size(mask)), source=0)

      j = 1
      do i=1, size(mask)
         if (mask(i) .ne. UNMARKED) then
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


   subroutine open_pointfile(filename, bdy_points, inside_pt, region_prefix)
      implicit none

      character(len=*) :: filename
      character(len=StrKIND), intent(out) :: region_prefix
      real(kind=RKIND), dimension(:,:), allocatable :: bdy_points
      real(kind=RKIND), dimension(2) :: inside_pt

      integer :: npts, i

      open(10, FILE=trim(filename))
      read(10, *) region_prefix
      read(10, *) npts
      allocate(bdy_points(2, npts))
      do i=1, npts
         read(10,*) bdy_points(1, i), bdy_points(2, i)
         bdy_points(:,i) = bdy_points(:,i) * PI / 180.0
      end do

      read(10,*) inside_pt(1), inside_pt(2)
      inside_pt = inside_pt * PI / 180.0

   end subroutine


!==================================================================================================
 integer function nearest_cell_path(target_lat, target_lon, start_cell, nCells, maxEdges, &
                               nEdgesOnCell, cellsOnCell, latCell, lonCell)
!==================================================================================================
 implicit none

 real (kind=RKIND), intent(in) :: target_lat, target_lon
 integer, intent(in) :: start_cell
 integer, intent(in) :: nCells, maxEdges
 integer, dimension(nCells), intent(in) :: nEdgesOnCell
 integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
 real (kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell

 integer :: i
 integer :: iCell
 integer :: current_cell
 real (kind=RKIND) :: current_distance, d
 real (kind=RKIND) :: nearest_distance

 nearest_cell_path = start_cell
 current_cell = -1

 do while (nearest_cell_path /= current_cell)
    current_cell = nearest_cell_path
    current_distance = sphere_distance(latCell(current_cell), lonCell(current_cell), target_lat, &
                                       target_lon, 1.0_RKIND)
    nearest_cell_path = current_cell
    nearest_distance = current_distance
    do i = 1, nEdgesOnCell(current_cell)
       iCell = cellsOnCell(i,current_cell)
       if (iCell <= nCells) then
          d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0_RKIND)
          if (d < nearest_distance) then
             nearest_cell_path = iCell
             nearest_distance = d
          end if
       end if
    end do
 end do

 end function nearest_cell_path

 function cross(u, v)
 implicit none

 real (kind=RKIND), dimension(3), intent(in) :: u, v
 real (kind=RKIND), dimension(3) :: cross

   cross(1) = u(2) * v(3) - v(2) * u(3)
   cross(2) = u(3) * v(1) - v(3) * u(1) 
   cross(3) = u(1) * v(2) - v(1) * u(2) 

 end function cross


 real(kind=RKIND) function dot(u, v)
 implicit none

 real (kind=RKIND), dimension(3), intent(in) :: u, v

 dot = u(1) * v(1) + u(2) * v(2) + u(3) * v(3)

 end function dot

 real(kind=RKIND) function mag(u)
 implicit none

 real (kind=RKIND), dimension(3), intent(in) :: u

 mag = u(1)**2 + u(2)**2 + u(3)**2
 mag = sqrt(mag)
 end function mag

end module utils
