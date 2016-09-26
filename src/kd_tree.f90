module kd_tree_mod
    use params
    implicit none
    type :: tree_node
        real(kind=RKIND), dimension(:), allocatable :: coords
        type(tree_node), pointer :: left_child => null()
        type(tree_node), pointer :: right_child => null()
        type(tree_node), pointer :: parent => null()
        integer :: id 
        contains 
        procedure :: print_node
        procedure :: nn_search
    end type

!    interface tree_node
!        module procedure treeconstructor
!    end interface
    
    type :: kd_tree
        type(tree_node), pointer :: root   
        
        contains
        procedure :: create_tree 
        procedure :: nearest_cell
    end type kd_tree

    contains

    integer function nearest_cell(this, q)
        class(kd_tree) :: this
        real(kind=RKIND), dimension(:), intent(in) :: q
        !integer, intent(out) :: nearest_cell
        integer :: k
        real(kind=RKIND) :: d

        d = huge(real(1.0, RKIND))
        k = size(q)
        call this%root%nn_search(1, k, nearest_cell, q, d)
    end function nearest_cell

    recursive subroutine nn_search(this, axs, k, id, q, dist)
        class(tree_node) :: this
        integer, intent(in) :: axs, k
        integer, intent(inout) :: id
        real(kind=RKIND), dimension(:), intent(in) :: q
        real(kind=RKIND), intent(inout) :: dist

        logical :: left_first
        real(kind=RKIND) :: temp
        left_first = q(axs) - this%coords(axs) <= 0.0

        if(left_first) then
            if(associated(this%left_child)) call this%left_child%nn_search(mod(axs,k)+1, k, id, q, dist)
            if(q(axs) + dist > this%coords(axs) .and. associated(this%right_child)) &
                call this%right_child%nn_search(mod(axs,k)+1, k, id, q, dist)
        else
            if(associated(this%right_child)) call this%right_child%nn_search(mod(axs,k)+1, k, id, q, dist)
            if(q(axs) - dist <= this%coords(axs) .and. associated(this%left_child)) &
                call this%left_child%nn_search(mod(axs,k)+1, k, id, q, dist)
        end if

        !temp = distance2(q, this%coords)
        temp = arc_distance(q, this%coords)
        if(temp < dist) then
            dist = temp
            id = this%id
        end if
    end subroutine nn_search
            
    function distance2(a, b)
        real(kind=RKIND), dimension(:), intent(in) :: a, b
        real(kind=RKIND) :: distance2
        
        integer :: k, i

        k = size(a)
        distance2 = 0.0
        do i=1,k
            distance2 = distance2 + (a(i) - b(i))**2
        end do
        distance2 = sqrt(distance2)
    end function distance2

    
    function arc_distance(a, b)
        real(kind=RKIND), dimension(:), intent(in) :: a, b
        real(kind=RKIND) :: arc_distance
    
        real(kind=RKIND) :: r2, angle
        integer :: k, i
        
        k = size(a)
        r2 = 0.0
        angle = 0.0
        do i=1, k
            r2 = r2 + a(i)**2
            angle = angle + a(i)*b(i)
        end do
        angle = angle / r2
        angle = acos(angle)
        r2=sqrt(r2)
        arc_distance = abs(angle*r2)
    end function arc_distance
        

    function keyCompare(a, b, l, k)
        real(kind=RKIND), dimension(:) :: a, b
        integer :: l, k
        integer :: i, j, c
        real(kind=RKIND) :: keyCompare
!        print *, "Enter keyCompare()"
        do i=l, k+l-1
            c = merge(i-k, i, i>k)
!            if (sum(a) == 0.0) print *, "ERROR: origin point as a"
!            if (sum(b) == 0.0) print *, "ERROR: origin point as b"
            keyCompare = a(c) - b(c)
!            print *, "in dim=", c, ",", a, "-", b, "=", keyCompare
            if(keyCompare .ne. 0.0) return
        end do
!        write (0,*) "ERROR: duplicate points in kd_tree: keyCompare()"
        return
    end function keyCompare 


    recursive subroutine mergeSort(points, ref, low, high, l, k)
        real(kind=RKIND), dimension(:,:) :: points
        integer, dimension(:), intent(inout) :: ref
        integer :: l, k, high, low

        integer :: i, j, m, n, med
        integer, dimension(:), allocatable :: temp
        
!        print *, "Enter mergeSort()"

        if(high > low) then

        n = high-low+1
        allocate(temp(low:high))
        med = n / 2 + low-1
        !if (n <= 2) med = low 
        call mergeSort(points, ref, low, med, l, k)
        call mergeSort(points, ref, med+1, high, l, k)

        temp(:) = ref(low:high)
        i = low
        j = med+1
!        print *, low, high, n, med
        do m=low, high
!            print *, "i=", i, "j=", j
            if (j <= high) then
                if (keyCompare(points(:,temp(i)), points(:,temp(j)), l, k) < 0 .and. i <= med) then
                    ref(m) = temp(i)
                    i = i+1
                end if
            else
                ref(m) = temp(j)
                j = j+1
            end if
        end do
        deallocate(temp)
        end if
!        print *, "Exit mergeSort()"
        
    end subroutine mergeSort


    recursive subroutine partition(points, ref, refout, low, high, n, l, k)
        real(kind=RKIND), dimension(:,:) :: points
        integer, dimension(:,:), intent(in) :: ref
        integer, dimension(:,:), allocatable :: parts
        integer, dimension(:), intent(inout) :: refout
        integer :: low, high, n, l, k, med, i, j, m, d, c
        
!        print*, "Enter partition(), n=", n
        med = (n + 1) / 2
        refout(low+med-1) = ref(med,l)
!        if (refout(low+med-1) == 0) print *, "Error: put a 0 in refout at ", low+med-1
!        print *, "med=", med, "n=", n, "k=", k
        if (n <= 1) return
        if (n == 2) then
            refout(low + med) = ref(med+1, l)
!            if (refout(low+med) == 0) print *, "Error(n==2): put a 0 in refout at ", low+med
            return
        end if
        if (n == 3) then
            refout(low) = ref(med-1, l)
            refout(high) = ref(med+1, l)
!            if (refout(low) == 0) print *, "Error(n==3): put a 0 in refout at ", low
!            if (refout(low+med-1) == 0) print *, "Error(n==3): put a 0 in refout at ", low+med-1
!            if (refout(high) == 0) print *, "Error(n==3): put a 0 in refout at ", high
            return
        end if
        allocate(parts(n, k))
        parts(:,l) = ref(:,l)
        do d=l+1, k+l
        c = merge(d-k, d, d>k)
        i = 1
        j = med+1
        do m=1, n
!            print *, "i=",i,"j=",j,"imax=", med-1, "jmax=", n
            if (ref(m, c) == refout(low+med-1)) cycle
            if (keyCompare(points(:, ref(m, c)), points(:, refout(low+med-1)), l, k) < 0 .and. i<med .or. j>n) then
                parts(i, c) = ref(m, c)
                i = i+1
            else
                parts(j, c) = ref(m, c)
                j = j+1
            end if
        end do
        end do

        call partition(points, parts(1:med-1, :), refout, low, low+med-2, med-1, mod(l, k)+1, k)
        call partition(points, parts(med+1:n, :), refout, low+med, high, n-med, mod(l, k)+1, k)

!        print *, "alloc:", allocated(parts), "med=", med, "n=", n, "k=", k
!        print *, parts
        !if(allocated(parts)) deallocate(parts)
!        print *, "Exit partition(), n=", n
    end subroutine partition

    subroutine sort(points, ids, n, k)
        real(kind=RKIND), dimension(:,:) :: points
        real(kind=RKIND), dimension(:,:), allocatable :: temp
        integer, dimension(:,:), allocatable :: ref
        integer, dimension(:), intent(inout) :: ids
        integer :: k, i, n, j
        
!        print *, "Enter sort()"
        allocate(temp(k, n))
        allocate(ref(n, k+1))
        do i=1,k
            ref(:, i) = (/(j, j=1, n)/)
            !if (i == 1) print *, ref(:,i)
            call mergeSort(points, ref(:,i), 1, n, i, k)
        end do
        call partition(points, ref(:,1:k), ref(:,k+1), 1, n, n, 1, k)
        temp = points
        do i=1,n
            points(:,i) = temp(:,ref(i,k+1))
        end do
        ids = ref(:,k+1)
!        print *, "Exit sort()"
    end subroutine sort
        

    subroutine create_tree(this, points)
        class(kd_tree) :: this
        real(kind=RKIND), dimension(:,:), intent(in) :: points
        real(kind=RKIND), dimension(:,:), allocatable :: sortpoints
        integer, dimension(:), allocatable :: ids
        integer :: n, k, i
        integer, dimension(2) :: dims

        dims = shape(points)
        k = dims(1)
        n = dims(2)
        allocate(sortpoints(k, n), source=points)
        allocate(ids(n))
        call sort(sortpoints, ids, n, k)            
!        print *, "Sorted Points:" 
!        do i=1, n
!            print *, ids(i), sortpoints(1,i), sortpoints(2,i), sortpoints(3,i)
!        end do
    
        this%root => build_tree(sortpoints, ids, n, k)
    end subroutine create_tree

    recursive function build_tree(points, ids, n, k) result(node)
        type(tree_node), pointer :: node
        real(kind=RKIND), dimension(:,:), intent(in) :: points
        integer, dimension(:), intent(in) :: ids
        integer :: median, n, k
    
        allocate(node)
        
        if(n == 1) then
            allocate(node%coords(k))
            node%coords(:)=points(:,1)
            node%id = ids(1)
            return
        end if    
        median = (n+1) / 2
        allocate(node%coords(k))
        node%coords(:) = points(:,median)
        node%id = ids(median)

        if(n > 2) node%left_child => build_tree(points(:,1:median-1), ids(1:median-1), median-1, k)
        if(n > 1) node%right_child => build_tree(points(:,median+1:n), ids(median+1:n), n-median, k)

    end function build_tree

    recursive subroutine print_node(this)
        class(tree_node) :: this
        if (associated(this%left_child)) call this%left_child%print_node()
        print *, new_line('A'), this%id, this%coords, new_line('A')
        if (associated(this%right_child)) call this%right_child%print_node()
    end subroutine print_node

end module kd_tree_mod

