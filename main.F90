program main

   use task_mod
   use hydro_data_mod
   USE river_lake_mod
   USE catchment_mod
   USE hillslope_mod
   USE information_mod
   USE output_mod

   implicit none

   INTEGER (kind = 4) :: narg
   character(len=256) :: hydro_dir, lake_dir, output_dir, casename, storage_type
   character(len=256) :: nlfile, rmfile
   real (kind=4) :: west,  east    ! -180 to 180
   real (kind=4) :: south, north   ! -90 to 90
   LOGICAL :: is_predefined_rivermouth
   
   real    (kind=4) :: catsize
   real    (kind=4) :: catsizemin   = 1.0
   real    (kind=4) :: lakecellsize = -1.e36
   integer (kind=4) :: nlev_max     = 10
   integer (kind=4) :: maxhunum
   integer (kind=4) :: maxnnb

   namelist /catexp/   &
      hydro_dir, lake_dir, output_dir, casename, storage_type, &
      catsize, lakecellsize, catsizemin, nlev_max, &
      west, east, south, north

   call task_init  ()

   if (p_is_master)  THEN
      narg = iargc ()

      CALL getarg (1, nlfile)
      open (10, status='OLD', file=nlfile, form="FORMATTED")
      read (10, nml=catexp)
      close(10)

      IF (narg >= 2) THEN
         is_predefined_rivermouth = .true.
         CALL getarg (2, rmfile)
      ELSE
         is_predefined_rivermouth = .false.
      ENDIF

      IF (lakecellsize <= 0) THEN
         lakecellsize = catsize / nlev_max
      ELSE
         lakecellsize = min(lakecellsize, catsize/nlev_max)
      ENDIF
   ENDIF

   call mpi_bcast (hydro_dir,    256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   call mpi_bcast (lake_dir,     256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   call mpi_bcast (output_dir,   256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   call mpi_bcast (casename,     256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   call mpi_bcast (storage_type, 256, MPI_CHARACTER, 0, p_comm_glb, p_err)
   call mpi_bcast (catsize,     1, MPI_REAL4, 0, p_comm_glb, p_err)
   call mpi_bcast (catsizemin,  1, MPI_REAL4, 0, p_comm_glb, p_err)
   call mpi_bcast (lakecellsize,1, MPI_REAL4, 0, p_comm_glb, p_err)
   call mpi_bcast (nlev_max, 1, MPI_INTEGER, 0, p_comm_glb, p_err)
   call mpi_bcast (west,  1, MPI_REAL4, 0, p_comm_glb, p_err)
   call mpi_bcast (east,  1, MPI_REAL4, 0, p_comm_glb, p_err)
   call mpi_bcast (south, 1, MPI_REAL4, 0, p_comm_glb, p_err)
   call mpi_bcast (north, 1, MPI_REAL4, 0, p_comm_glb, p_err)
   call mpi_bcast (is_predefined_rivermouth, 1, MPI_LOGICAL, 0, p_comm_glb, p_err)
   
   ! Step 0: Initializing data and work.
   if (p_is_master) write(*,'(A)') 'Step 0 : Data and Work specialization ...'
   call get_region (hydro_dir, west, east, north, south)
   call divide_data_and_work (count(bkid==0))

   if (p_is_master) write(*,'(/A)') 'Step 0 : Loading input data ...'
   call readin_blocks (hydro_dir, lake_dir)

   call init_area ()
   CALL init_window ()
   
   ! Step 1: Finding all rivers and lake outlets.
   IF (is_predefined_rivermouth) THEN
      call get_river_lake (catsize, catsizemin, rmfile)
   ELSE
      call get_river_lake (catsize, catsizemin)
   ENDIF

   ! Step 2: Dividing the region into catchments.
   call get_catchment (catsize)
      
   ! Step 3: Dividing catchment into hillslopes and hydrounits.
   call get_hillslope_hydrounits (catsize, lakecellsize, nlev_max, maxhunum)

   ! Step 4: Get river and hillslope information.
   call get_information (maxhunum, maxnnb)

   ! Step 5: Writing out results.
   CALL output_result (output_dir, casename, storage_type, maxhunum, maxnnb)
     
   call free_memory ()

   call task_final ()

end program main
