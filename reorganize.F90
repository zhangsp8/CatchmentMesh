PROGRAM reorgnize 

   USE task_mod
   USE hydro_data_mod
   USE river_lake_mod
   USE catchment_mod
   USE hillslope_mod
   USE neighbour_mod
   USE output_mod
   USE ncio_serial

   IMPLICIT NONE

   character(len=256) :: nlfile, rmfile, filename
   real (kind=4) :: west,  east    ! -180 to 180
   real (kind=4) :: south, north   ! -90 to 90
   logical :: end_of_data
   integer :: iblk, jblk

   real    (kind=4) :: catsize
   real    (kind=4) :: catsizemin   = 1.0
   real    (kind=4) :: lakecellsize = -1.e36
   integer (kind=4) :: nlev_max     = 10
   integer (kind=4) :: nhrumax
   integer (kind=4) :: nnbmax
   integer (kind=4) :: ntotalall

   namelist /catexp/   &
      hydro_dir, lake_dir, output_dir, casename, storage_type, &
      catsize, lakecellsize, catsizemin, nlev_max, &
      west, east, south, north

   CALL task_init  ()

   IF (p_is_master)  THEN
      CALL getarg (1, nlfile)
      open (10, status='OLD', file=nlfile, form="FORMATTED")
      read (10, nml=catexp)
      close(10)

      IF (lakecellsize <= 0) THEN
         lakecellsize = catsize / nlev_max
      ELSE
         lakecellsize = min(lakecellsize, catsize/nlev_max)
      ENDIF
   ENDIF

   CALL mpi_bcast (hydro_dir,    256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (lake_dir,     256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (output_dir,   256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (casename,     256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (storage_type, 256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (catsize,        1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (catsizemin,     1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (lakecellsize,   1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (nlev_max,       1,   MPI_INTEGER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (west,           1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (east,           1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (south,          1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (north,          1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   
   ! Step 0: Initializing data and work.
   IF (p_is_master) write(*,'(A)') 'Step 0 : Assign Data and Work processors ...'

   CALL get_region (west, east, north, south)
   CALL init_area ()

   thisinfo => allinfo
   thisinfo%icatdsp = 0
   thisinfo%next => null()

   ntotalall = 0

   nhrumax = 0
   
   DO WHILE (.true.)

      CALL get_next_rivermouth (catsize, catsizemin, end_of_data)
      IF (end_of_data) THEN
         EXIT
      ENDIF

      IF (p_is_master) THEN

         iblk = thisinfo%ithisblk
         jblk = thisinfo%jthisblk
         
         CALL get_filename (trim(output_dir) // '/' // trim(casename), iblk, jblk, filename)
         
         write(*,*) 'Read block results from file ', trim(filename)

         CALL ncio_read_serial (filename, 'bsn_index', thisinfo%bsn_index)
         CALL ncio_read_serial (filename, 'lake_id'  , thisinfo%lake_id  )
         CALL ncio_read_serial (filename, 'bsn_elva' , thisinfo%bsn_elva )
         CALL ncio_read_serial (filename, 'bsn_downstream', thisinfo%bsn_downstream)
         CALL ncio_read_serial (filename, 'bsn_nswe' , thisinfo%bsn_nswe )

         CALL ncio_read_serial (filename, 'riv_len', thisinfo%riv_len)
         CALL ncio_read_serial (filename, 'riv_elv', thisinfo%riv_elv)
         
         CALL ncio_read_serial (filename, 'bsn_num_hru', thisinfo%bsn_num_hru)
         
         CALL ncio_read_serial (filename, 'hru_indx' , thisinfo%hru_indx)
         CALL ncio_read_serial (filename, 'hru_area' , thisinfo%hru_area)
         CALL ncio_read_serial (filename, 'hru_hand' , thisinfo%hru_hand)
         CALL ncio_read_serial (filename, 'hru_elva' , thisinfo%hru_elva)
         CALL ncio_read_serial (filename, 'hru_next' , thisinfo%hru_next)
         CALL ncio_read_serial (filename, 'hru_plen' , thisinfo%hru_plen)
         CALL ncio_read_serial (filename, 'hru_lfac' , thisinfo%hru_lfac)

         thisinfo%ntotalcat = size(thisinfo%bsn_index)
         nhrumax = max(nhrumax, size(thisinfo%hru_indx,1))

         write(*,*) 'Number of catchments on ', iblk, jblk, ' is ', thisinfo%ntotalcat
         write(*,*) 'Number of hydrounits on ', iblk, jblk, ' is ', size(thisinfo%hru_indx,1)

      ENDIF

   ENDDO

   ! Step 4: Get basin neighbour information.
   CALL get_basin_neighbour (ntotalall, nnbmax)

   ! Step 5: Writing out results.
   CALL output_result (nhrumax, nnbmax)
     
   CALL free_memory ()

END PROGRAM
