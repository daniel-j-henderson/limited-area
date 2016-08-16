module minheap_mod
   use params
   implicit none
   type :: heap_node
      integer :: cellID, parent_cell
      real (kind=RKIND) :: priority
      contains
         procedure :: setNodeEqual
         generic :: assignment(=) => setNodeEqual
   end type heap_node

   type :: min_heap
      integer :: max_cap, curr_cap, size
      integer, dimension(:), allocatable :: index_array
      type(heap_node), dimension(:), allocatable :: h
      contains
         procedure :: is_empty
         procedure :: create_heap
         procedure :: delete_heap
         procedure :: decrease_priority
         procedure :: extract_min
         procedure :: insert
         procedure :: reheapify_up
         procedure :: reheapify_down
   end type min_heap


   contains

   subroutine setNodeEqual(this, a)
      class(heap_node), intent(out) :: this
      type(heap_node), intent(in) :: a
      this%priority = a%priority
      this%cellID = a%cellID
      this%parent_cell = a%parent_cell
   end subroutine setNodeEqual

   function is_empty(this)
      class(min_heap) :: this
      logical :: is_empty
      if (this%size == 0) then
         is_empty = .true.
      else
         is_empty = .false.
      end if
   end function

   subroutine create_heap(this, init_cap, max_cap, id, p)
      implicit none
      class(min_heap) :: this
      integer :: init_cap, max_cap, id
      real (kind=RKIND) :: p
      
      this%curr_cap = init_cap
      this%max_cap = max_cap
      allocate(this%index_array(max_cap), source=0)
      allocate(this%h(this%curr_cap))
      this%h(1)%priority = p
      this%h(1)%cellID = id
      this%h(1)%parent_cell = id
      this%index_array(id) = 1
      this%size = 1
      
   end subroutine create_heap

   subroutine delete_heap(this)
      implicit none
      class(min_heap) :: this
      
      if (allocated(this%h)) deallocate(this%h)
      if (allocated(this%index_array)) deallocate(this%index_array)
      this%size = 0
   end subroutine delete_heap
   
   
   subroutine decrease_priority(this, id, p)
      implicit none
      class(min_heap) :: this
      integer :: id
      real (kind=RKIND) :: p
   
      integer :: my_index, parent_index, parent_id
      real (kind=RKIND) my_priority, parent_priority 
   
      if (this%index_array(id) == 0) then
         write (0,*) "ERROR: Trying to decrease priority of an element which is not present"
         return
      end if
      
      ! change priority value
      this%h(this%index_array(id))%priority = p
      
      ! re-heapify up
      call this%reheapify_up(this%index_array(id))

   end subroutine decrease_priority

   function extract_min(this)
      implicit none
      class(min_heap) :: this
      integer :: extract_min

      integer :: my_id

      extract_min = this%h(1)%cellID

      if (this%size == 1) then
         this%size = 0
         return
      end if
    
      ! put the last item in the root position and decrease the size by 1
      my_id = this%h(this%size)%cellID !my_id id cellID of last element
      this%index_array(my_id) = 1
      this%h(1) = this%h(this%size)
      this%h(1)%parent_cell = my_id
      this%size = this%size - 1  

      !re-heapify down

      call this%reheapify_down(1)
   end function

   subroutine insert(this, my_id, priority)
      implicit none
      class(min_heap) :: this
      integer :: my_id
      real (kind=RKIND) :: priority

      type(heap_node), dimension(:), allocatable :: node_temp
      integer, dimension(:), allocatable :: index_temp
      
      
      if (.not. this%size < this%curr_cap) then
         if (this%curr_cap < this%max_cap) then
            this%curr_cap = min(this%curr_cap * 2, this%max_cap)
            allocate(node_temp(size(this%h)), source=this%h)
            deallocate(this%h)  !, this%index_array)
            allocate(this%h(this%curr_cap))  !, this%index_array(this%curr_cap))
            this%h(1:this%size) = node_temp(1:this%size)
         else
            write (0,*) "Cannot insert new item, heap is at capacity"
         end if
      end if


      this%size = this%size+1
      this%index_array(my_id) = this%size
      this%h(this%size)%priority = priority
      this%h(this%size)%cellID = my_id
      this%h(this%size)%parent_cell = this%h(this%size / 2)%cellID

      call this%reheapify_up(this%size)

   end subroutine insert

   subroutine reheapify_up(this, from_index)
      implicit none
      class(min_heap) :: this
      integer :: from_index

      integer :: my_index, parent_index, my_id, parent_id
      real (kind=RKIND) :: my_priority, parent_priority

      my_index = from_index
      my_id = this%h(my_index)%cellID
      parent_id = this%h(my_index)%parent_cell
      parent_index = this%index_array(parent_id)
      my_priority = this%h(my_index)%priority
      parent_priority = this%h(parent_index)%priority
      do while (my_priority < parent_priority)
         ! give my higher index (worse priority) to my parent
         ! take my parent's better index
         this%index_array(my_id) = parent_index
         this%index_array(parent_id) = my_index

         ! at my parent's location insert my data
         this%h(parent_index)%priority = my_priority
         this%h(parent_index)%cellID = my_id
      
         ! at my location insert my parent's data
         this%h(my_index)%priority = parent_priority
         this%h(my_index)%cellID = parent_id
         this%h(my_index)%parent_cell = my_id
         ! now i have swapped places with my parent, who is a worse priority than me (higher real value = lower priority)

         ! now prepare for new compare
         my_index = this%index_array(my_id)
         parent_id = this%h(my_index)%parent_cell
         parent_index = this%index_array(parent_id)
         my_priority = this%h(my_index)%priority
         parent_priority = this%h(parent_index)%priority

         !TODO - be sure to make sure everything works when a node overtakes the root
      end do
   end subroutine reheapify_up

   subroutine reheapify_down(this, from_index)
      implicit none
      class(min_heap) :: this
      integer :: from_index

      integer :: my_index, lc_index, rc_index, swap_index, my_id, swap_id
      real (kind=RKIND) :: my_priority, lc_priority, rc_priority, swap_priority

      my_index = from_index
      lc_index = 2*my_index
      rc_index = lc_index+1
      my_id = this%h(my_index)%cellID
      my_priority = this%h(my_index)%priority
      lc_priority = merge(this%h(lc_index)%priority, huge(1.0_RKIND), lc_index <= this%size)
      rc_priority = merge(this%h(rc_index)%priority, huge(1.0_RKIND), rc_index <= this%size)

      do while(my_priority > lc_priority .or. my_priority > rc_priority)
         if (lc_priority < rc_priority) then
            ! swap with left child
            swap_index = lc_index
         else
            ! swap with right child
            swap_index = rc_index
         end if

         swap_id = this%h(swap_index)%cellID
         swap_priority = this%h(swap_index)%priority

         this%index_array(swap_id) = my_index
         this%index_array(my_id) = swap_index

         this%h(swap_index)%priority = my_priority
         this%h(swap_index)%cellID = my_id
         this%h(swap_index)%parent_cell = swap_id

         this%h(my_index)%priority = swap_priority
         this%h(my_index)%cellID = swap_id

         my_index = this%index_array(my_id)
         lc_index = 2*my_index
         rc_index = lc_index+1
         my_priority = this%h(my_index)%priority
         lc_priority = merge(this%h(lc_index)%priority, huge(1.0_RKIND), lc_index <= this%size)
         rc_priority = merge(this%h(rc_index)%priority, huge(1.0_RKIND), rc_index <= this%size)
      end do
   end subroutine reheapify_down

end module minheap_mod
