PROGRAM get_rivermouth

   USE task_mod
   USE hydro_data_mod
   IMPLICIT NONE
   
   character(len=256) :: hydro_dir, output_dir, casename, storage_type
   character(len=256) :: nlfile
   real (kind=4) :: west,  east    ! -180 to 180
   real (kind=4) :: south, north   ! -90 to 90
   real    (kind=4) :: catsize
   integer (kind=4) :: nlev_max

   integer (kind=4) :: i, j
   integer (kind=4) :: i_up, j_up, i_dn, j_dn
   integer (kind=4) :: iblk, jblk
   
   integer (kind=4) :: iproc, ndone
   real    (kind=4) :: upa

   namelist /catexp/   &
      hydro_dir, output_dir, casename, storage_type, &
      catsize, nlev_max, west, east, south, north

   CALL task_init  ()

   IF (p_is_master)  THEN
      CALL getarg (1, nlfile)
      open (10, status='OLD', file=nlfile, form="FORMATTED")
      read (10, nml=catexp)
      close(10)
   ENDIF

   CALL mpi_bcast (hydro_dir,    256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   CALL mpi_bcast (output_dir,   256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   CALL mpi_bcast (casename,     256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   CALL mpi_bcast (storage_type, 256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   CALL mpi_bcast (catsize,  1, MPI_REAL4,   0, p_comm_glb, p_err)
   CALL mpi_bcast (nlev_max, 1, MPI_INTEGER, 0, p_comm_glb, p_err)
   CALL mpi_bcast (west,  1, MPI_REAL4, 0, p_comm_glb, p_err)
   CALL mpi_bcast (east,  1, MPI_REAL4, 0, p_comm_glb, p_err)
   CALL mpi_bcast (south, 1, MPI_REAL4, 0, p_comm_glb, p_err)
   CALL mpi_bcast (north, 1, MPI_REAL4, 0, p_comm_glb, p_err)
   
   ! Step 0: Initializing data and work.
   IF (p_is_master) write(*,*) 'Step 0: Data and Work specialization ...'
   CALL get_region (hydro_dir, west, east, north, south)
   CALL divide_data_and_work (count(bkid==0))

   IF (p_is_master) write(*,'(/A)') 'Step 0: Loading input data ...'
   CALL readin_blocks (hydro_dir)

   CALL init_area ()

   CALL mpi_barrier (p_comm_glb)

   IF (p_is_data) THEN

      DO jblk = 1, mblock
         DO iblk = 1, nblock
            IF (bkid(iblk,jblk) == p_iam_glb) THEN
               DO j = 1, mbox
                  DO i = 1, nbox
                     IF (blks(iblk,jblk)%icat(i,j) == 0) THEN
                        ! river mouth and inland depression.
                        IF ((blks(iblk,jblk)%dir(i,j) == 0) .or. (blks(iblk,jblk)%dir(i,j) == -1)) THEN
                           i_up = blks(iblk,jblk)%idsp + i
                           j_up = blks(iblk,jblk)%jdsp + j
                           CALL mpi_send (p_iam_glb, 1, MPI_INTEGER, 0, 0, p_comm_glb, p_err) 
                           CALL mpi_send (i_up, 1, MPI_INTEGER, 0, 1, p_comm_glb, p_err) 
                           CALL mpi_send (j_up, 1, MPI_INTEGER, 0, 1, p_comm_glb, p_err) 
                           CALL mpi_send (blks(iblk,jblk)%upa(i,j), 1, MPI_REAL4, 0, 1, p_comm_glb, p_err) 
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
                              
      CALL mpi_send (-1, 1, MPI_INTEGER, 0, 0, p_comm_glb, p_err) 

   ELSEIF (p_is_master) THEN

      ndone = 0
      DO WHILE (.true.)
         CALL mpi_recv (iproc, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0, p_comm_glb, p_stat, p_err)

         IF (iproc > 0) THEN
            CALL mpi_recv (i_up,  1, MPI_INTEGER, iproc, 1, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (j_up,  1, MPI_INTEGER, iproc, 1, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (upa,   1, MPI_REAL4,   iproc, 1, p_comm_glb, p_stat, p_err)
            write (*,'(I7,I7,E20.6)') i_up, j_up, upa
         ELSE
            ndone = ndone + 1
         ENDIF

         IF (ndone == p_ndata) EXIT
      ENDDO

   ENDIF
   
   CALL free_memory ()

   CALL task_final ()

END PROGRAM get_rivermouth
