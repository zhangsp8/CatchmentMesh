PROGRAM get_rivermouth

   use task_mod
   use hydro_data_mod
   implicit none
   
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

   call task_init  ()

   if (p_is_master)  THEN
      CALL getarg (1, nlfile)
      open (10, status='OLD', file=nlfile, form="FORMATTED")
      read (10, nml=catexp)
      close(10)
   ENDIF

   call mpi_bcast (hydro_dir,    256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   call mpi_bcast (output_dir,   256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   call mpi_bcast (casename,     256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   call mpi_bcast (storage_type, 256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   call mpi_bcast (catsize,  1, MPI_REAL4,   0, p_comm_glb, p_err)
   call mpi_bcast (nlev_max, 1, MPI_INTEGER, 0, p_comm_glb, p_err)
   call mpi_bcast (west,  1, MPI_REAL4, 0, p_comm_glb, p_err)
   call mpi_bcast (east,  1, MPI_REAL4, 0, p_comm_glb, p_err)
   call mpi_bcast (south, 1, MPI_REAL4, 0, p_comm_glb, p_err)
   call mpi_bcast (north, 1, MPI_REAL4, 0, p_comm_glb, p_err)
   
   ! Step 0: Initializing data and work.
   if (p_is_master) write(*,*) 'Step 0: Data and Work specialization ...'
   call get_region (hydro_dir, west, east, north, south)
   call divide_data_and_work (count(bkid==0))

   if (p_is_master) write(*,'(/A)') 'Step 0: Loading input data ...'
   call readin_blocks (hydro_dir)

   call init_area ()

   CALL mpi_barrier (p_comm_glb)

   if (p_is_data) then

      do jblk = 1, mblock
         do iblk = 1, nblock
            if (bkid(iblk,jblk) == p_iam_glb) then
               do j = 1, mbox
                  do i = 1, nbox
                     if (blks(iblk,jblk)%icat(i,j) == 0) then
                        ! river mouth and inland depression.
                        IF ((blks(iblk,jblk)%dir(i,j) == 0) .or. (blks(iblk,jblk)%dir(i,j) == -1)) THEN
                           i_up = blks(iblk,jblk)%idsp + i
                           j_up = blks(iblk,jblk)%jdsp + j
                           call mpi_send (p_iam_glb, 1, MPI_INTEGER, 0, 0, p_comm_glb, p_err) 
                           call mpi_send (i_up, 1, MPI_INTEGER, 0, 1, p_comm_glb, p_err) 
                           call mpi_send (j_up, 1, MPI_INTEGER, 0, 1, p_comm_glb, p_err) 
                           call mpi_send (blks(iblk,jblk)%upa(i,j), 1, MPI_REAL4, 0, 1, p_comm_glb, p_err) 
                        end if
                     end if
                  end do
               end do
            end if
         end do
      end do
                              
      call mpi_send (-1, 1, MPI_INTEGER, 0, 0, p_comm_glb, p_err) 

   elseif (p_is_master) THEN

      ndone = 0
      DO WHILE (.true.)
         call mpi_recv (iproc, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0, p_comm_glb, p_stat, p_err)

         IF (iproc > 0) THEN
            call mpi_recv (i_up,  1, MPI_INTEGER, iproc, 1, p_comm_glb, p_stat, p_err)
            call mpi_recv (j_up,  1, MPI_INTEGER, iproc, 1, p_comm_glb, p_stat, p_err)
            call mpi_recv (upa,   1, MPI_REAL4,   iproc, 1, p_comm_glb, p_stat, p_err)
            write (*,'(I7,I7,E20.6)') i_up, j_up, upa
         ELSE
            ndone = ndone + 1
         ENDIF

         IF (ndone == p_ndata) exit
      ENDDO

   end if
   
   call free_memory ()

   call task_final ()

end PROGRAM get_rivermouth
