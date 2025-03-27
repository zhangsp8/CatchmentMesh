PROGRAM main

   USE task_mod
   USE hydro_data_mod
   USE river_lake_mod
   USE catchment_mod
   USE hillslope_mod
   USE neighbour_mod
   USE output_mod

   IMPLICIT NONE

   integer (kind = 4) :: narg
   character(len=256) :: deldir, nlfile, rmfile, filename
   real (kind=4) :: west  = -180.
   real (kind=4) :: east  =  180.
   real (kind=4) :: south = -90.
   real (kind=4) :: north =  90.
   logical :: has_predefined_rivermouth, end_of_data

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
      west, east, south, north, include_all_upstream

   CALL task_init  ()

   IF (p_is_master)  THEN
      narg = iargc ()

      CALL getarg (1, nlfile)
      open (10, status='OLD', file=nlfile, form="FORMATTED")
      read (10, nml=catexp)
      close(10)

      IF (narg >= 2) THEN
         has_predefined_rivermouth = .true.
         CALL getarg (2, rmfile)
      ELSE
         has_predefined_rivermouth = .false.
      ENDIF

      IF (lakecellsize <= 0) THEN
         lakecellsize = catsize / nlev_max
      ELSE
         lakecellsize = min(lakecellsize, catsize/nlev_max)
      ENDIF
   ENDIF

   CALL mpi_bcast (hydro_dir,           256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (lake_dir,            256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (output_dir,          256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (casename,            256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (storage_type,        256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (catsize,               1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (catsizemin,            1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (lakecellsize,          1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (nlev_max,              1,   MPI_INTEGER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (west,                  1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (east,                  1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (south,                 1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (north,                 1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (include_all_upstream,  1,   MPI_LOGICAL, p_master_address, p_comm_glb, p_err)

   CALL mpi_bcast (has_predefined_rivermouth, 1, MPI_LOGICAL, p_master_address, p_comm_glb, p_err)
   
   IF (has_predefined_rivermouth) THEN
      storage_type = 'one'
   ENDIF

   ! Step 0: Initializing data and work.
   IF (p_is_master) write(*,'(A)') 'Step 0 : Assign Data and Work processors ...'

   IF (p_is_master .and. (trim(storage_type) == 'block')) THEN
      deldir = 'rm -rf ' // trim(output_dir) // '/' // trim(casename)
      CALL system(trim(deldir))
   ENDIF

   CALL get_region (west, east, north, south)
   CALL init_area ()

   thisinfo => allinfo
   thisinfo%icatdsp = 0
   thisinfo%next => null()

   ntotalall = 0

   nhrumax = 1
   
   DO WHILE (.true.)

      ! Step 1: Finding all rivers and lake outlets.
      IF (has_predefined_rivermouth) THEN

         CALL get_river_lake (catsize, catsizemin, rivermouthfile = rmfile)

      ELSE

         CALL get_river_lake (catsize, catsizemin, end_of_data = end_of_data)
      
         IF (end_of_data) THEN
            EXIT
         ENDIF
      ENDIF

      ! Step 2: Dividing the region into catchments.
      CALL get_catchment (catsize)

      ! Step 3: Dividing catchment into hillslopes and hydrounits.
      CALL get_hillslope_hydrounits (catsize, lakecellsize, nlev_max, nhrumax)

      IF (p_is_master) THEN
         ntotalall = ntotalall + thisinfo%ntotalcat
      ENDIF
      CALL mpi_bcast (ntotalall, 1, MPI_INTEGER, p_master_address, p_comm_glb, p_err)

      IF (p_is_master .and. (trim(storage_type) == 'block')) THEN
         CALL output_block_info (nhrumax)
         CALL flush_blocks (output = .true.)
      ENDIF
         
      IF (has_predefined_rivermouth) THEN
         EXIT
      ENDIF

   ENDDO

   ! Step 4: Get basin neighbour information.
   CALL get_basin_neighbour (ntotalall, nnbmax)

   ! Step 5: Writing out results.
   CALL output_result (nhrumax, nnbmax)
     
   CALL free_memory ()

END PROGRAM main
