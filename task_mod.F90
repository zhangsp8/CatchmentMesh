MODULE task_mod

   ! USE mpi
   IMPLICIT NONE
   
   include 'mpif.h'

   logical :: p_is_master    
   logical :: p_is_work
   integer :: p_master_address
   integer :: p_nwork

   ! Global communicator information
   integer :: p_comm_glb
   integer :: p_iam_glb    
   integer :: p_nglb

   integer :: p_stat(MPI_STATUS_SIZE)
   integer :: p_err

CONTAINS

   !-------------------------------------
   SUBROUTINE task_init ()

      IMPLICIT NONE
      integer :: errorcode
   
      CALL mpi_init (p_err) 

      p_comm_glb = MPI_COMM_WORLD
      CALL mpi_comm_rank (p_comm_glb, p_iam_glb, p_err)  
      CALL mpi_comm_size (p_comm_glb, p_nglb, p_err) 

      p_master_address = 0

      p_is_master = (p_iam_glb == p_master_address)
      p_is_work   = (p_iam_glb /= p_master_address)

      p_nwork = p_nglb-1

   END SUBROUTINE task_init

END MODULE task_mod
