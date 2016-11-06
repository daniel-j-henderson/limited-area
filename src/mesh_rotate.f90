 module mesh_rotate
 use params

   real(kind=RKIND) :: pii = PI

 contains

   subroutine rotate(xCell, yCell, zCell, original_lat, original_lon, new_lat, new_lon)

      implicit none
      real (kind=RKIND), pointer, dimension(:), intent(inout) :: xCell, yCell, zCell
      real (kind=RKIND) :: original_lat, original_lon, new_lat, new_lon

      real (kind=RKIND) :: xNew, yNew, zNew, thetaLat, thetaLon, x0LongitudeAtEquator, y0LongitudeAtEquator, z0LongitudeAtEquator, uCrossProduct, vCrossProduct, wCrossProduct
      integer :: i

      thetaLat = new_lat- original_lat
      thetaLon = new_lon- original_lon

      ! create the unit vector <x0LongitudeAtEquator,
      ! y0LongitudeAtEquator, z0LongitudeAtEquator>
      call convert_lx(x0LongitudeAtEquator, y0LongitudeAtEquator, z0LongitudeAtEquator, 1.0_RKIND, 0.0_RKIND, original_lon)

      ! create the unit vector <xNew, yNew, zNew>
      call convert_lx(xNew, yNew, zNew, 1.0_RKIND, new_lat, new_lon)

      ! create the unit vector <uCrossProduct, vCrossProduct,
      ! wCrossProduct> by using a right-angle cross-product of two unit
      ! vectors in the perpendicular plane
      call cross_product(x0LongitudeAtEquator, y0LongitudeAtEquator, z0LongitudeAtEquator, &
                         0.0_RKIND, 0.0_RKIND, 1.0_RKIND, &
                         uCrossProduct, vCrossProduct, wCrossProduct)

      do i=1,size(xCell)
         call executeRotation(xCell(i), yCell(i), zCell(i), thetaLat, thetaLon, uCrossProduct, vCrossProduct, wCrossProduct, xNew, yNew, zNew)
      end do

   end subroutine rotate 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !  Note that the since <u,v,w> must be a unit vector (where u**2 +
   !  v**2 + w**2 = 1),
   !  <uCrossProduct, vCrossProduct, wCrossProduct> and <xNew,yNew,zNew>
   !  must be unit vectors as well
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine executeRotation(x, y, z, thetaLat, thetaLon, uCrossProduct, vCrossProduct, wCrossProduct, xNew, yNew, zNew)

      implicit none

         real (kind=RKIND), intent(inout) :: x, y, z
         real (kind=RKIND), intent(in) :: thetaLat, thetaLon
         real (kind=RKIND), intent(in) :: uCrossProduct, vCrossProduct, wCrossProduct
         real (kind=RKIND), intent(in) :: xNew, yNew, zNew
         real (kind=RKIND) ::  u, v, w

         ! latitude rotation (rotate around cross product of xyz
         ! corresponding to original point's longitude at the equator
         ! and the z axis unit vector)
         u = uCrossProduct
         v = vCrossProduct
         w = wCrossProduct
         call rotate_about_vector(x, y, z, thetaLat, 0.0_RKIND, 0.0_RKIND, 0.0_RKIND, u, v, w)

         ! longitude rotation (rotate around z axis)
         u = 0.0
         v = 0.0
         w = 1.0
         call rotate_about_vector(x, y, z, thetaLon, 0.0_RKIND, 0.0_RKIND, 0.0_RKIND, u, v, w)


   end subroutine executeRotation


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CONVERT_XL
   !
   ! Convert (x, y, z) to a (lat, lon) location on a sphere with
   !    radius sqrt(x^2 + y^2 + z^2).
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine convert_xl(x, y, z, lat, lon)
   
      implicit none
   
      real (kind=RKIND), intent(in) :: x, y, z
      real (kind=RKIND), intent(out) :: lat, lon
   
      real (kind=RKIND) :: dl, clat
      real (kind=RKIND) :: eps
      parameter (eps=1.e-10)
   
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

   end subroutine convert_xl


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CONVERT_LX
   !
   ! Convert (lat,lon) to an (x, y, z) location on a sphere with
   ! specified radius.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine convert_lx(x, y, z, radius, lat, lon)

      implicit none

      real (kind=RKIND), intent(in) :: radius, lat, lon
      real (kind=RKIND), intent(out) :: x, y, z

      z = radius * sin(lat)
      x = radius * cos(lon) * cos(lat)
      y = radius * sin(lon) * cos(lat)

   end subroutine convert_lx


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION SPHERE_ANGLE
   !
   ! Computes the angle between arcs AB and AC, given points A, B, and C
   ! Equation numbers w.r.t.
   ! http://mathworld.wolfram.com/SphericalTrigonometry.html
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real function sphere_angle(ax, ay, az, bx, by, bz, cx, cy, cz)

      implicit none

      real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz

      real (kind=RKIND) :: a, b, c          ! Side lengths of spherical triangle ABC

      real (kind=RKIND) :: ABx, ABy, ABz    ! The components of the vector AB
      real (kind=RKIND) :: mAB              ! The magnitude of AB
      real (kind=RKIND) :: ACx, ACy, ACz    ! The components of the vector AC
      real (kind=RKIND) :: mAC              ! The magnitude of AC

      real (kind=RKIND) :: Dx               ! The i-components of the cross product AB x AC
      real (kind=RKIND) :: Dy               ! The j-components of the cross product AB x AC
      real (kind=RKIND) :: Dz               ! The k-components of the cross product AB x AC

      real (kind=RKIND) :: s                ! Semiperimeter of the triangle
      real (kind=RKIND) :: sin_angle

      a = acos(max(min(bx*cx + by*cy + bz*cz,1.0_RKIND),-1.0_RKIND))      ! Eqn. (3)
      b = acos(max(min(ax*cx + ay*cy + az*cz,1.0_RKIND),-1.0_RKIND))      ! Eqn. (2)
      c = acos(max(min(ax*bx + ay*by + az*bz,1.0_RKIND),-1.0_RKIND))      ! Eqn. (1)

      ABx = bx - ax
      ABy = by - ay
      ABz = bz - az

      ACx = cx - ax
      ACy = cy - ay
      ACz = cz - az

      Dx =   (ABy * ACz) - (ABz * ACy)
      Dy = -((ABx * ACz) - (ABz * ACx))
      Dz =   (ABx * ACy) - (ABy * ACx)

      s = 0.5*(a + b + c)
!      sin_angle = sqrt((sin(s-b)*sin(s-c))/(sin(b)*sin(c)))   ! Eqn.
!      (28)
      sin_angle = sqrt(min(1.0_RKIND,max(0.0_RKIND,(sin(s-b)*sin(s-c))/(sin(b)*sin(c)))))   ! Eqn. (28)

      if ((Dx*ax + Dy*ay + Dz*az) >= 0.0) then
         sphere_angle =  2.0 * asin(max(min(sin_angle,1.0_RKIND),-1.0_RKIND))
      else
         sphere_angle = -2.0 * asin(max(min(sin_angle,1.0_RKIND),-1.0_RKIND))
      end if

   end function sphere_angle
   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ROTATE_ABOUT_VECTOR
   !
   ! Rotates the point (x,y,z) through an angle theta about the line
   ! through (a,b,c) 
   ! with direction vector <u,v,w>. Note that the uvw must describe a
   ! unit vector (where u**2 + v**2 + w**2 = 1) 
   !
   ! Reference:
   ! http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine rotate_about_vector(x, y, z, theta, a, b, c, u, v, w)

      implicit none

      real (kind=RKIND), intent(inout) :: x, y, z
      real (kind=RKIND), intent(in) :: theta, a, b, c, u, v, w
      real (kind=RKIND) :: xp, yp, zp

      real (kind=RKIND) :: vw2, uw2, uv2
      real (kind=RKIND) :: m

      vw2 = v**2.0 + w**2.0
      uw2 = u**2.0 + w**2.0
      uv2 = u**2.0 + v**2.0
      m = sqrt(u**2.0 + v**2.0 + w**2.0)

      xp = (a*vw2 + u*(-b*v-c*w+u*x+v*y+w*z) + ((x-a)*vw2+u*(b*v+c*w-v*y-w*z))*cos(theta) + m*(-c*v+b*w-w*y+v*z)*sin(theta))/m**2.0
      yp = (b*uw2 + v*(-a*u-c*w+u*x+v*y+w*z) + ((y-b)*uw2+v*(a*u+c*w-u*x-w*z))*cos(theta) + m*( c*u-a*w+w*x-u*z)*sin(theta))/m**2.0
      zp = (c*uv2 + w*(-a*u-b*v+u*x+v*y+w*z) + ((z-c)*uv2+w*(a*u+b*v-u*x-v*y))*cos(theta) + m*(-b*u+a*v-v*x+u*y)*sin(theta))/m**2.0

      ! alternate calculation
      !xp = (a*vw2 - u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(theta)) +
      !x*cos(theta) + (-c*v+b*w-w*y+v*z)*sin(theta)
      !yp = (b*uw2 - v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(theta)) +
      !y*cos(theta) + ( c*u-a*w+w*x-u*z)*sin(theta)
      !zp = (c*uv2 - w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(theta)) +
      !z*cos(theta) + (-b*u+a*v-v*x+u*y)*sin(theta)

      x = xp
      y = yp
      z = zp

   end subroutine rotate_about_vector


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! SUBROUTINE CROSS_PRODUCT
   !
   ! Computes C = A x B
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine cross_product(Ax, Ay, Az, &
                            Bx, By, Bz, &
                            Cx, Cy, Cz)

      real (kind=RKIND), intent(in)  :: Ax, Ay, Az
      real (kind=RKIND), intent(in)  :: Bx, By, Bz
      real (kind=RKIND), intent(out) :: Cx, Cy, Cz

      Cx = (Ay * Bz) - (Az * By)
      Cy = (Az * Bx) - (Ax * Bz)
      Cz = (Ax * By) - (Ay * Bx)

   end subroutine cross_product                                 

end module 
