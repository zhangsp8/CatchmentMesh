MODULE task_mod

   ! USE mpi
   IMPLICIT NONE
   
   include 'mpif.h'

   logical :: p_is_master    
   logical :: p_is_data
   logical :: p_is_work

   integer :: p_master_address
   integer :: p_data_address

   ! Global communicator information
   integer :: p_comm_glb
   integer :: p_iam_glb    
   integer :: p_nglb

   ! Input/output and data processors 
   integer :: p_comm_data
   integer :: p_iam_data
   integer :: p_ndata     

   ! Processors carrying out the computing work
   integer :: p_comm_work
   integer :: p_iam_work
   integer :: p_nwork     

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
      p_data_address   = 1

      p_is_master = (p_iam_glb == p_master_address)
      p_is_data   = (p_iam_glb == p_data_address)
      p_is_work   = (.not. p_is_master) .and. (.not. p_is_data)

      p_ndata = 1
      p_nwork = p_nglb-2


      IF (p_is_data .or. p_is_master) THEN
         CALL mpi_comm_split (p_comm_glb, 1, p_iam_glb, p_comm_data, p_err)
         CALL mpi_comm_rank  (p_comm_data, p_iam_data, p_err)  
      ELSE
         CALL mpi_comm_split (p_comm_glb, MPI_UNDEFINED, p_iam_glb, p_comm_data, p_err)
      ENDIF

      IF (p_is_work .or. p_is_master) THEN
         CALL mpi_comm_split (p_comm_glb, 1, p_iam_glb, p_comm_work, p_err)
         
         CALL mpi_comm_size (p_comm_work, p_nwork,   p_err) 
         CALL mpi_comm_rank (p_comm_work, p_iam_work, p_err)  
         p_nwork = p_nwork - 1
      ELSE
         CALL mpi_comm_split (p_comm_glb, MPI_UNDEFINED, p_iam_glb, p_comm_work, p_err)
      ENDIF
      
   END SUBROUTINE task_init

END MODULE task_mod
