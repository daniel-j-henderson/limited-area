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




   subroutine create_boundary_from_region(bdyMaskCell, regionType, rparams, &
                                          latCell, lonCell, radius, cellsOnCell, nEdgesOnCell)

      implicit none
   
      integer, dimension(:), pointer, intent(inout) :: bdyMaskCell
      integer, dimension(:), pointer, intent(in) :: nEdgesOnCell
      integer, dimension(:,:), pointer, intent(in) :: cellsOnCell
      integer, intent(in) :: regionType
      real(kind=RKIND), dimension(:), intent(in) :: rparams
      real(kind=RKIND), dimension(:), pointer, intent(in) :: latCell, lonCell 
      real(kind=RKIND), intent(in) :: radius

      integer :: i, j, nCells
      real(kind=RKIND) x, y, min_axis_angle, lat_center, lon_center, lat_a, lon_a, lat_b, lon_b, lat_f1, lon_f1, lat_f2, lon_f2, phi, temp, com_dist
      real(kind=RKIND), dimension(3) :: r, r2, a_pt, b_pt, c_pt, f1, f2

      bdyMaskCell = UNMARKED
      nCells = size(bdyMaskCell)
      select case (regionType)
      case (RCIRCULAR)
         ! rparams = (/lat, lon, circle_region_radius/)
         if (rparams(3) > PI * radius) then
            write (0,*) "Error: The radius you provided for your circular region is greater than the radius of the MPAS mesh * PI. Please provide a new radius."
            stop
         end if
         lat_center = rparams(1) * PI / 180.0
         lon_center = rparams(2) * PI / 180.0
         do i=1, nCells
            if (sphere_distance(lat_center, lon_center, latCell(i), lonCell(i), radius) <= rparams(3)) &
               bdyMaskCell(i) = INSIDE
         end do

      case (RELLIPTICAL)
         ! rparams = (/lat_center, lon_center, lat_a, lon_a, minor_axis_length_degrees/)   
         lat_center = rparams(1) * PI / 180.0
         lon_center = rparams(2) * PI / 180.0
         lat_a = rparams(3) * PI / 180.0
         lon_a = rparams(4) * PI / 180.0
         min_axis_angle = rparams(5) * PI / 180.0
         !lat_b = rparams(5) * PI / 180.0
         !lon_b = rparams(6) * PI / 180.0
         call con_lx(lat_a, lon_a, radius, a_pt(1), a_pt(2), a_pt(3))
         call con_lx(lat_b, lon_b, radius, b_pt(1), b_pt(2), b_pt(3))
         call con_lx(lat_center, lon_center, radius, c_pt(1), c_pt(2), c_pt(3))

         ! rotate the center point about the unit normal vector by the
         ! appropriate amount to find the foci
         r = cross(a_pt, c_pt)
         r = r/mag(r)

         ! b is calculated so that you need not provide it, only its angle
         r2 = cross(a_pt, r)
         r2 = r/mag(r)
         b_pt = rot(c_pt, r2, min_axis_angle)
         call con_xl(b_pt(1), b_pt(2), b_pt(3), lat_b, lon_b)

         phi = sqrt(sphere_distance(lat_a, lon_a, lat_center, lon_center, radius)**2 - &
                    (min_axis_angle * radius)**2)
         phi = phi / radius
         f1 = rot(c_pt, r, -phi)
         f2 = rot(c_pt, r, phi)
         
         call con_xl(f1(1), f1(2), f1(3), lat_f1, lon_f1)
         call con_xl(f2(1), f2(2), f2(3), lat_f2, lon_f2)

         ! distance from one focus to any point on the ellipse to the other
         ! focus, we compare with that to see which points lie inside
         com_dist = sphere_distance(lat_f1, lon_f1, lat_a, lon_a, radius) + &
                    sphere_distance(lat_f2, lon_f2, lat_a, lon_a, radius)
         do i=1, nCells
            x = sphere_distance(latCell(i), lonCell(i), lat_f1, lon_f1, radius)
            y = sphere_distance(latCell(i), lonCell(i), lat_f2, lon_f2, radius)
            if ((x + y) <= com_dist) then
                  bdyMaskCell(i) = INSIDE
            end if
         end do
      case default
      end select

      call mark_neighbors_of_type(INSIDE, BOUNDARY1, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY1, BOUNDARY2, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY2, BOUNDARY3, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY3, BOUNDARY4, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY4, BOUNDARY5, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY5, BOUNDARY6, bdyMaskCell, cellsOnCell, nEdgesOnCell)
      call mark_neighbors_of_type(BOUNDARY6, BOUNDARY7, bdyMaskCell, cellsOnCell, nEdgesOnCell)


   end subroutine create_boundary_from_region


    subroutine create_boundary_custom(nCells, radius, bdy_pts, inside_pt, bdyMaskCell, nEdgesOnCell, cellsOnCell, latCell, lonCell)


      implicit none

      integer :: iCell, source_cell, v, i, j, k, npts, target_cell, inside_cell
      integer, intent(in) :: nCells
      real (kind=RKIND), dimension(:), pointer, intent(in) :: latCell, lonCell
      real (kind=RKIND), dimension(:,:), intent(in) :: bdy_pts
      real (kind=RKIND), dimension(2) :: inside_pt
      integer, dimension(:), pointer, intent(out) :: bdyMaskCell
      integer, dimension(:), pointer, intent(in) :: nEdgesOnCell
      integer, dimension(:,:), pointer, intent(in) :: cellsOnCell
      integer, dimension(:), allocatable :: boundary_cells
      real (kind=RKIND) :: mindist, dist, angle, minangle 
      real (kind=RKIND), intent(in) :: radius
      real (kind=RKIND), dimension(3) :: pt, pta, ptb, v1, v2, v3

      real(kind=RKIND) :: temp
      


      ! Interpolate the boundary points provided onto the MPAS mesh
      npts = size(bdy_pts(1,:))
      allocate(boundary_cells(npts))
      do i=1, size(bdy_pts(1,:))
         boundary_cells(i) = nearest_cell_path(bdy_pts(1,i), bdy_pts(2,i), &
                  merge(boundary_cells(i-1), 1, i > 1), nCells, 10, nEdgesOnCell, cellsOnCell, latCell, lonCell)
      end do
      inside_cell = nearest_cell_path(inside_pt(1), inside_pt(2), &
                  boundary_cells(1), nCells, 10, nEdgesOnCell, cellsOnCell, latCell, lonCell)
         
      bdyMaskCell = UNMARKED

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

   end subroutine create_boundary_custom

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


   subroutine open_pointfile(filename, bdy_points, inside_pt, region_prefix, region_type, rparams)
      implicit none

      character(len=*) :: filename
      character(len=StrKIND), intent(out) :: region_prefix
      character(len=StrKIND) :: rtype 
      integer, intent(out) :: region_type
      character(len=1000) :: str
      real(kind=RKIND), dimension(:,:), allocatable :: bdy_points
      real(kind=RKIND), dimension(2) :: inside_pt
      real(kind=RKIND), dimension(:), allocatable, intent(out) :: rparams

      integer :: npts, i, j, n

      open(10, FILE=trim(filename))
      read(10, *) region_prefix
      read(10, *) str
      if (verify(trim(str), '0123456789') > 0) then
         read(str, *) rtype
         if (trim(rtype) == 'circle') then
            region_type = RCIRCULAR
            n = 3
         else if (trim(rtype) == 'ellipse') then
            region_type = RELLIPTICAL
            n = 5
         end if
         allocate(rparams(n))
         read(10, *) rparams
      else 
         region_type = RCUSTOM
         read(str, *) npts
         allocate(bdy_points(2, npts))
         do i=1, npts
            read(10,*) bdy_points(1, i), bdy_points(2, i)
            bdy_points(:,i) = bdy_points(:,i) * PI / 180.0
         end do

         read(10,*) inside_pt(1), inside_pt(2)
         inside_pt = inside_pt * PI / 180.0
      end if

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


 function rot(v, k, phi)
   implicit none
   real(kind=RKIND), dimension(3) :: rot
   real(kind=RKIND), dimension(3), intent(in) :: v, k
   real(kind=RKIND) :: phi

   rot = v*cos(phi) + cross(k, v)*sin(phi) + v*dot(k, v)*(1-cos(phi))

 end function rot


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CONVERT_XL
   !
   ! Convert (x, y, z) to a (lat, lon) location on a sphere with
   !    radius sqrt(x^2 + y^2 + z^2).
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine con_xl(x, y, z, lat, lon)
   
      implicit none
   
      real (kind=RKIND), intent(in) :: x, y, z
      real (kind=RKIND), intent(out) :: lat, lon
   
      real (kind=RKIND) :: dl
      real (kind=RKIND) :: clat, eps
      parameter (eps=1.e-10)
      
      integer :: i,j
      real(kind=RKIND) :: pii = PI

   
      dl = sqrt(x*x + y*y + z*z)
      lat = asin(z/dl)
   
   !  check for being close to either pole
         if (abs(x) > eps) then
            if (abs(y) > eps) then
               lon = atan(abs(y/x))
      
               if ((x <= 0.) .and. (y >= 0.)) then
                  lon = pii-lon
               else if ((x <= 0.) .and. (y < 0.)) then
                  lon = lon+pii
               else if ((x >= 0.) .and. (y <= 0.)) then
                  lon = 2*pii-lon
               end if
      
            else ! we're either on longitude 0 or 180
      
               if (x > 0) then
                  lon = 0.
               else
                  lon = pii
               end if
      
            end if
      
         else if (abs(y) > eps) then
      
            if (y > 0) then
               lon = pii/2.
            else
               lon = 3.*pii/2.
            end if
      
         else  ! we are at a pole
      
            lon = 0.
      
         end if
   end subroutine con_xl


end module utils
