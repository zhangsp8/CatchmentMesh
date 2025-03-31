PROGRAM main

   USE task_mod
   USE hydro_data_mod
   USE river_lake_mod
   USE catchment_mod
   USE hillslope_mod
   USE neighbour_mod
   USE output_mod

   IMPLICIT NONE

   real (kind=4) :: def_west  = -180.
   real (kind=4) :: def_east  =  180.
   real (kind=4) :: def_south = -90.
   real (kind=4) :: def_north =  90.

   real    (kind=4) :: def_catsize
   real    (kind=4) :: def_catsizemin   = 1.0
   real    (kind=4) :: def_lakecellsize = -1.e36
   integer (kind=4) :: def_num_profile  = 10
   integer (kind=4) :: def_nlev_max     = 10

   integer (kind = 4) :: narg
   character(len=256) :: deldir, nlfile, rmfile, filename
   logical :: has_predefined_rivermouth, end_of_data

   integer (kind=4) :: nhrumax
   integer (kind=4) :: nnbmax
   integer (kind=4) :: ntotalall

   namelist /catexp/   &
      def_hydro_dir,    def_lake_dir,    def_output_dir,   def_casename,   &
      def_storage_type, def_catsize,     def_lakecellsize, def_catsizemin, &
      def_west,         def_east,        def_south,        def_north,      &
      def_nlev_max,     def_num_profile, def_include_all_upstream

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

      IF (def_lakecellsize <= 0) THEN
         def_lakecellsize = def_catsize / def_nlev_max
      ELSE
         def_lakecellsize = min(def_lakecellsize, def_catsize/def_nlev_max)
      ENDIF
   ENDIF

   CALL mpi_bcast (def_hydro_dir,           256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_lake_dir,            256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_output_dir,          256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_casename,            256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_storage_type,        256, MPI_CHARACTER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_catsize,               1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_catsizemin,            1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_lakecellsize,          1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_west,                  1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_east,                  1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_south,                 1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_north,                 1,     MPI_REAL4, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_nlev_max,              1,   MPI_INTEGER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_num_profile,           1,   MPI_INTEGER, p_master_address, p_comm_glb, p_err)
   CALL mpi_bcast (def_include_all_upstream,  1,   MPI_LOGICAL, p_master_address, p_comm_glb, p_err)

   CALL mpi_bcast (has_predefined_rivermouth, 1, MPI_LOGICAL, p_master_address, p_comm_glb, p_err)

   IF (has_predefined_rivermouth) THEN
      def_storage_type = 'one'
   ENDIF

   ! Step 0: Initializing data and work.
   IF (p_is_master) write(*,'(A)') 'Step 0 : Assign Data and Work processors ...'

   IF (p_is_master .and. (trim(def_storage_type) == 'block')) THEN
      deldir = 'rm -rf ' // trim(def_output_dir) // '/' // trim(def_casename)
      CALL system(trim(deldir))
   ENDIF

   CALL get_region (def_west, def_east, def_north, def_south)
   CALL init_area ()

   thisinfo => allinfo
   thisinfo%icatdsp = 0
   thisinfo%next => null()

   thisinfo%ithisblk = -1
   thisinfo%jthisblk = -1

   ntotalall = 0

   nhrumax = 1

   DO WHILE (.true.)

      ! Step 1: Finding all rivers and lake outlets.
      IF (has_predefined_rivermouth) THEN

         CALL get_river_lake (def_catsize, def_catsizemin, rivermouthfile = rmfile)

      ELSE

         CALL get_river_lake (def_catsize, def_catsizemin, end_of_data = end_of_data)

         IF (end_of_data) THEN
            EXIT
         ENDIF
      ENDIF

      ! Step 2: Dividing the region into catchments.
      CALL get_catchment (def_catsize)

      ! Step 3: Dividing catchment into hillslopes and hydrounits.
      CALL get_hillslope_hydrounits (def_catsize, def_lakecellsize, def_nlev_max, nhrumax)

      IF (p_is_master) THEN
         ntotalall = ntotalall + thisinfo%ntotalcat
      ENDIF
      CALL mpi_bcast (ntotalall, 1, MPI_INTEGER, p_master_address, p_comm_glb, p_err)

      IF (p_is_master .and. (trim(def_storage_type) == 'block')) THEN
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
