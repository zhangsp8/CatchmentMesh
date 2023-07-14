module task_mod

   implicit none
   
   include 'mpif.h'

   ! Global communicator information
   integer :: p_comm_glb
   integer :: p_iam_glb    
   integer :: p_nglb
   logical :: p_is_master    

   ! Input/output and data processors 
   integer :: p_comm_data
   integer :: p_iam_data
   integer :: p_ndata     
   logical :: p_is_data

   ! Processors carrying out the computing work
   integer :: p_comm_work
   integer :: p_iam_work
   integer :: p_nwork     
   logical :: p_is_work

   ! Processors sharing the same memory
   integer :: p_comm_node
   integer :: p_iam_node
   integer :: p_nnode

   logical :: has_master_node
   integer :: p_ndata_node

   integer :: p_stat(MPI_STATUS_SIZE)
   integer :: p_err

   integer, allocatable :: p_data_list (:)
   integer, allocatable :: p_inode_all (:)

contains

   !-------------------------------------
   subroutine task_init ()

      implicit none
   
      call mpi_init (p_err) 

      p_comm_glb = MPI_COMM_WORLD
      call mpi_comm_rank (p_comm_glb, p_iam_glb, p_err)  
      call mpi_comm_size (p_comm_glb, p_nglb, p_err) 

      p_is_master = (p_iam_glb == 0)

      call mpi_comm_split_type (p_comm_glb, MPI_COMM_TYPE_SHARED, p_iam_glb, &
         MPI_INFO_NULL, p_comm_node, p_err)
      call mpi_comm_rank (p_comm_node, p_iam_node, p_err)  
      call mpi_comm_size (p_comm_node, p_nnode, p_err)

      if (p_nnode < 3) then 
         write(*,*) 'Please use at least 3 processors on each node.'
         call mpi_abort ()
      end if

   end subroutine task_init

   !-------------------------------------
   subroutine divide_data_and_work (ndata)

      implicit none
   
      integer, intent(in) :: ndata
   
      logical, allocatable :: p_is_all (:)
      integer :: inode, nnode, ithis, idata, iproc

      p_ndata = min(ndata, (p_nglb-1)/2)
   
      nnode = 0
      ithis = p_nglb
      do while (ithis > 0)
         if ((p_iam_node == 0) .and. (p_iam_glb < ithis)) then
            ithis = p_iam_glb
         else
            ithis = 0
         end if

         call mpi_allreduce (MPI_IN_PLACE, ithis, 1, MPI_INTEGER, MPI_MAX, p_comm_glb, p_err)

         nnode = nnode + 1
         if (p_iam_glb == ithis) inode = nnode
      end do
   
      call mpi_bcast (inode, 1, MPI_INTEGER, 0, p_comm_node, p_err)

      p_ndata_node = p_ndata / nnode
      if (inode <= mod(p_ndata,nnode)) then
         p_ndata_node = p_ndata_node + 1
      end if

      has_master_node = p_is_master
      call mpi_allreduce (MPI_IN_PLACE, has_master_node, 1, MPI_LOGICAL, MPI_LOR, p_comm_node, p_err)

      if (has_master_node) then
         p_is_data = (.not. p_is_master) .and. (p_iam_node <= p_ndata_node)
      else
         p_is_data = (p_iam_node < p_ndata_node)
      end if

      p_is_work = (.not. p_is_master) .and. (.not. p_is_data)

      if (p_is_data .or. p_is_master) then
         call mpi_comm_split (p_comm_glb, 1, p_iam_glb, p_comm_data, p_err)
         call mpi_comm_rank  (p_comm_data, p_iam_data, p_err)  
      else
         call mpi_comm_split (p_comm_glb, MPI_UNDEFINED, p_iam_glb, p_comm_data, p_err)
      end if

      if (p_is_work .or. p_is_master) then
         call mpi_comm_split (p_comm_glb, 1, p_iam_glb, p_comm_work, p_err)
         
         call mpi_comm_size (p_comm_work, p_nwork,   p_err) 
         call mpi_comm_rank (p_comm_work, p_iam_work, p_err)  
         p_nwork = p_nwork - 1
      else
         call mpi_comm_split (p_comm_glb, MPI_UNDEFINED, p_iam_glb, p_comm_work, p_err)
      end if
      
      call mpi_bcast (p_nwork, 1, MPI_INTEGER, 0, p_comm_glb, p_err)
   
      allocate (p_is_all (0:p_nglb-1))
      call mpi_allgather (p_is_data, 1, MPI_LOGICAL, p_is_all, 1, MPI_LOGICAL, &
         p_comm_glb, p_err)

      p_ndata = count(p_is_all)
      allocate (p_data_list (p_ndata))
      idata = 0
      do iproc = 0, p_nglb-1
         if (p_is_all(iproc)) then
            idata = idata + 1
            p_data_list(idata) = iproc
         end if
      end do
      
      allocate (p_inode_all (0:p_nglb-1))
      call mpi_allgather (inode, 1, MPI_INTEGER, p_inode_all, 1, MPI_INTEGER, p_comm_glb, p_err)

      if (p_is_master) then

         write(*,100) p_iam_glb, p_inode_all(0)
         100 format('Master is','(',I4,','I4,')') 

         write(*,*) 'Data in on:'
         do iproc = 0, p_nglb-1
            if (p_is_all(iproc)) then
               write(*,200) iproc, p_inode_all(iproc)
               200 format('(',I4,','I4,')',',',$)
            end if
         end do
         write(*,*)

         write(*,*) 'Work in on:'
         do iproc = 0, p_nglb-1
            if ((.not. p_is_all(iproc)) .and. (iproc /= p_iam_glb)) then
               write(*,300) iproc, p_inode_all(iproc)
               300 format('(',I4,','I4,')',',',$)
            end if
         end do
         write(*,*)
      end if
      
      deallocate (p_is_all)
      
   end subroutine divide_data_and_work 

   !------------------------------------
   subroutine task_final ()

      implicit none

      if (allocated(p_data_list))   deallocate (p_data_list)
      if (allocated(p_inode_all))   deallocate (p_inode_all)

      call mpi_finalize (p_err)

   end subroutine task_final

end module task_mod
