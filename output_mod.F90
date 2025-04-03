MODULE output_mod

CONTAINS

   SUBROUTINE output_block_info (nhrumax, nstep)

      USE task_mod
      USE hydro_data_mod
      USE ncio_serial

      IMPLICIT NONE

      integer (kind=4), intent(in) :: nhrumax
      integer (kind=4), intent(in) :: nstep

      integer :: iblk, jblk
      character(len=256) :: mkdirname, filename
      logical :: fexists

      IF (p_is_master) THEN

         iblk = thisinfo%ithisblk
         jblk = thisinfo%jthisblk

         mkdirname = 'mkdir -p ' // trim(def_output_dir) // '/' // trim(def_casename)
         CALL system(trim(mkdirname))

         CALL get_filename (trim(def_output_dir) // '/' // trim(def_casename), iblk, jblk, filename)

         write(*,*) 'Write block results to file ', trim(filename)

         inquire (file=trim(filename), exist=fexists)

         IF (.not. fexists) THEN
            CALL ncio_create_file (filename)
            CALL ncio_define_dimension (filename, 'lat',  nbox)
            CALL ncio_define_dimension (filename, 'lon', mbox)
            CALL ncio_write_serial (filename, 'lat', blks(iblk,jblk)%lat, 'lat')
            CALL ncio_write_serial (filename, 'lon', blks(iblk,jblk)%lon, 'lon')
            CALL ncio_write_serial (filename, 'elva', blks(iblk,jblk)%elv,  &
               'lat', 'lon', compress = 1)
         ENDIF

         IF (thisinfo%ntotalcat > 0) THEN

            CALL ncio_define_dimension (filename, 'basin', thisinfo%ntotalcat)
            CALL ncio_define_dimension (filename, 'nswe' , 4)
            CALL ncio_define_dimension (filename, 'step' , nstep)

            CALL ncio_write_serial (filename, 'bsn_index', thisinfo%bsn_index, 'basin', compress = 1)
            CALL ncio_write_serial (filename, 'lake_id'  , thisinfo%lake_id  , 'basin', compress = 1)
            CALL ncio_write_serial (filename, 'bsn_elva' , thisinfo%bsn_elva,  'basin', compress = 1)
            CALL ncio_write_serial (filename, 'bsn_downstream' , thisinfo%bsn_downstream,  'basin', compress = 1)
            CALL ncio_write_serial (filename, 'bsn_nswe' , thisinfo%bsn_nswe,  'nswe', 'basin', compress = 1)

            CALL ncio_define_dimension (filename, 'river', thisinfo%nrivseg)
            CALL ncio_write_serial (filename, 'riv_len', thisinfo%riv_len, 'river', compress = 1)
            CALL ncio_write_serial (filename, 'riv_elv', thisinfo%riv_elv, 'river', compress = 1)

            CALL ncio_write_serial (filename, 'bsn_num_hru', thisinfo%bsn_num_hru, 'basin', compress = 1)

            CALL ncio_define_dimension (filename, 'hydrounit' , nhrumax)
            CALL ncio_write_serial (filename, 'hru_indx' , thisinfo%hru_indx,  'hydrounit', 'basin', compress = 1)
            CALL ncio_write_serial (filename, 'hru_area' , thisinfo%hru_area,  'hydrounit', 'basin', compress = 1)
            CALL ncio_write_serial (filename, 'hru_hand' , thisinfo%hru_hand,  'hydrounit', 'basin', compress = 1)
            CALL ncio_write_serial (filename, 'hru_elva' , thisinfo%hru_elva,  'hydrounit', 'basin', compress = 1)
            CALL ncio_write_serial (filename, 'hru_next' , thisinfo%hru_next,  'hydrounit', 'basin', compress = 1)
            CALL ncio_write_serial (filename, 'hru_plen' , thisinfo%hru_plen,  'hydrounit', 'basin', compress = 1)
            CALL ncio_write_serial (filename, 'hru_lfac' , thisinfo%hru_lfac,  'hydrounit', 'basin', compress = 1)

            CALL ncio_write_serial (filename, 'hru_fdstep' , thisinfo%hru_fdstep,  &
               'step', 'hydrounit', 'basin', compress = 1)

         ENDIF

      ENDIF

   END SUBROUTINE output_block_info

   SUBROUTINE output_result (nhrumax, nnbmax, nstep)

      USE task_mod
      USE hydro_data_mod
      USE ncio_serial

      IMPLICIT NONE

      integer (kind=4), intent(in) :: nhrumax
      integer (kind=4), intent(in) :: nnbmax
      integer (kind=4), intent(in) :: nstep

      integer (kind=4) :: np, mp
      real    (kind=8), allocatable :: longitude(:)
      real    (kind=8), allocatable :: latitude (:)
      integer (kind=4), allocatable :: catch (:,:)
      real    (kind=4), allocatable :: hnd   (:,:)
      real    (kind=4), allocatable :: elv   (:,:)
      integer (kind=4), allocatable :: hunit (:,:)

      type(info_typ) :: outinfo

      character(len=256) :: filename, mkdirname
      integer :: nthis, dsp, nhru, iunit, iblk, jblk, iloc, ic, ndim1
      logical :: found
      real (kind = 4), parameter :: spval = -1.e36

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) write(*,'(/A)') 'Step 5: Output results ...'


      IF (p_is_master) THEN

         IF (trim(def_storage_type) == 'one') THEN

            ! shrink regions
            inorth = allinfo%bsn_nswe(1,1)
            isouth = allinfo%bsn_nswe(2,1)
            jwest  = allinfo%bsn_nswe(3,1)
            jeast  = allinfo%bsn_nswe(4,1)

            thisinfo => allinfo
            DO WHILE (associated(thisinfo))

               DO ic = 1, thisinfo%ntotalcat
                  inorth = min     (inorth, thisinfo%bsn_nswe(1,ic))
                  isouth = max     (isouth, thisinfo%bsn_nswe(2,ic))
                  jwest  = min_west(jwest , thisinfo%bsn_nswe(3,ic))
                  jeast  = max_east(jeast , thisinfo%bsn_nswe(4,ic))
               ENDDO

               thisinfo => thisinfo%next
            ENDDO

            write(*,*) 'region after shrink (nswe): ', inorth, isouth, jwest, jeast

            filename = trim(def_output_dir) // '/' // trim(def_casename) // '.nc'

            CALL ncio_create_file  (filename)
            CALL ncio_write_serial (filename, 'inorth', inorth)
            CALL ncio_write_serial (filename, 'isouth', isouth)
            CALL ncio_write_serial (filename, 'jwest', jwest)
            CALL ncio_write_serial (filename, 'jeast', jeast)

            np = isouth - inorth + 1
            mp = jeast  - jwest  + 1
            IF (mp < 0) mp = mp + mglb

            allocate (latitude (np))
            allocate (longitude(mp))
            allocate (catch (np,mp))
            allocate (hnd   (np,mp))
            allocate (elv   (np,mp))
            allocate (hunit (np,mp))

            CALL aggregate_data (inorth, isouth, jwest, jeast, &
               np, mp, longitude = longitude, latitude = latitude, &
               icat = catch, hnd = hnd, elv = elv, hunit = hunit)

            CALL ncio_define_dimension (filename, 'lat', np)
            CALL ncio_define_dimension (filename, 'lon', mp)
            CALL ncio_write_serial (filename, 'lat', latitude,  'lat')
            CALL ncio_write_serial (filename, 'lon', longitude, 'lon')

            CALL ncio_put_attr_str (filename, 'lat', 'long_name', 'latitude')
            CALL ncio_put_attr_str (filename, 'lat', 'units', 'degrees_north')
            CALL ncio_put_attr_str (filename, 'lon', 'long_name', 'longitude')
            CALL ncio_put_attr_str (filename, 'lon', 'units', 'degrees_east')


            CALL ncio_write_serial (filename, 'icatchment2d', catch, 'lat', 'lon', compress = 1)
            CALL ncio_write_serial (filename, 'hand', hnd, 'lat', 'lon', compress = 1)
            CALL ncio_write_serial (filename, 'elva', elv, 'lat', 'lon', compress = 1)
            CALL ncio_write_serial (filename, 'ihydrounit2d', hunit, 'lat', 'lon', compress = 1)

            deallocate (latitude )
            deallocate (longitude)
            deallocate (catch)
            deallocate (hnd  )
            deallocate (elv  )
            deallocate (hunit)

         ENDIF

         outinfo%ntotalcat = 0
         thisinfo => allinfo
         DO WHILE (associated(thisinfo))
            outinfo%ntotalcat = outinfo%ntotalcat + thisinfo%ntotalcat
            thisinfo => thisinfo%next
         ENDDO

         filename = trim(def_output_dir) // '/' // trim(def_casename) // '.nc'
         IF (trim(def_storage_type) == 'block') THEN
            CALL ncio_create_file (filename)
         ENDIF

         CALL ncio_define_dimension (filename, 'catchment', outinfo%ntotalcat)
         CALL ncio_define_dimension (filename, 'hydrounit', nhrumax)
         CALL ncio_define_dimension (filename, 'neighbour', nnbmax)
         CALL ncio_define_dimension (filename, 'step', nstep)

         ! ----- output : number of hydro units in a basin -----
         allocate(outinfo%bsn_num_hru (outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) outinfo%bsn_num_hru(dsp+1:dsp+nthis) = thisinfo%bsn_num_hru
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'basin_numhru', outinfo%bsn_num_hru, &
            'catchment', compress = 1)

         ! ----- output : hydro unit index -----
         allocate(outinfo%hru_indx (nhrumax, outinfo%ntotalcat)); outinfo%hru_indx(:,:) = -1
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) nhru  = size(thisinfo%hru_indx,1)
            IF (nthis > 0) outinfo%hru_indx(1:nhru,dsp+1:dsp+nthis) = thisinfo%hru_indx
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_index', outinfo%hru_indx, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit area -----
         allocate(outinfo%hru_area (nhrumax, outinfo%ntotalcat)); outinfo%hru_area(:,:) = 0.
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) nhru  = size(thisinfo%hru_area,1)
            IF (nthis > 0) outinfo%hru_area(1:nhru,dsp+1:dsp+nthis) = thisinfo%hru_area
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_area', outinfo%hru_area, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit hand -----
         allocate(outinfo%hru_hand (nhrumax, outinfo%ntotalcat)); outinfo%hru_hand(:,:) = 0.
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) nhru  = size(thisinfo%hru_hand,1)
            IF (nthis > 0) outinfo%hru_hand(1:nhru,dsp+1:dsp+nthis) = thisinfo%hru_hand
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_hand', outinfo%hru_hand, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit elevation -----
         allocate(outinfo%hru_elva (nhrumax, outinfo%ntotalcat)); outinfo%hru_elva(:,:) = 0.
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) nhru  = size(thisinfo%hru_elva,1)
            IF (nthis > 0) outinfo%hru_elva(1:nhru,dsp+1:dsp+nthis) = thisinfo%hru_elva
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_elva', outinfo%hru_elva, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit downstream -----
         allocate(outinfo%hru_next (nhrumax, outinfo%ntotalcat)); outinfo%hru_next(:,:) = -1
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) nhru  = size(thisinfo%hru_next,1)
            IF (nthis > 0) outinfo%hru_next(1:nhru,dsp+1:dsp+nthis) = thisinfo%hru_next
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_downstream', outinfo%hru_next, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit path length -----
         allocate(outinfo%hru_plen (nhrumax, outinfo%ntotalcat)); outinfo%hru_plen(:,:) = 0.
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) nhru  = size(thisinfo%hru_plen,1)
            IF (nthis > 0) outinfo%hru_plen(1:nhru,dsp+1:dsp+nthis) = thisinfo%hru_plen
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_pathlen', outinfo%hru_plen, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit drianage contour length -----
         allocate(outinfo%hru_lfac (nhrumax, outinfo%ntotalcat)); outinfo%hru_lfac(:,:) = 0.
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) nhru  = size(thisinfo%hru_lfac,1)
            IF (nthis > 0) outinfo%hru_lfac(1:nhru,dsp+1:dsp+nthis) = thisinfo%hru_lfac
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_facelen', outinfo%hru_lfac, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit flood area profile -----
         allocate(outinfo%hru_fdstep (nstep, nhrumax, outinfo%ntotalcat)); outinfo%hru_fdstep(:,:,:) = 0.
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) nhru  = size(thisinfo%hru_fdstep,2)
            IF (nthis > 0) outinfo%hru_fdstep(:,1:nhru,dsp+1:dsp+nthis) = thisinfo%hru_fdstep
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_flood_step', outinfo%hru_fdstep, &
            'step', 'hydrounit', 'catchment', compress = 1)

         ! ----- output : river length -----
         allocate(outinfo%riv_len (outinfo%ntotalcat));   outinfo%riv_len(:) = spval
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%nrivseg
            IF (nthis > 0) outinfo%riv_len(dsp+1:dsp+nthis) = thisinfo%riv_len
            dsp = dsp + thisinfo%ntotalcat
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'river_length', outinfo%riv_len, &
            'catchment', compress = 1)
         CALL ncio_put_attr_real4 (filename, 'river_length', 'missing_value', spval)

         ! ----- output : river elevation -----
         allocate(outinfo%riv_elv (outinfo%ntotalcat));   outinfo%riv_elv(:) = spval
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%nrivseg
            IF (nthis > 0) outinfo%riv_elv(dsp+1:dsp+nthis) = thisinfo%riv_elv
            dsp = dsp + thisinfo%ntotalcat
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'river_elevation', outinfo%riv_elv, &
            'catchment', compress = 1)
         CALL ncio_put_attr_real4 (filename, 'river_elevation', 'missing_value', spval)

         ! ----- output : downstream basin -----
         allocate(outinfo%bsn_downstream (outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) outinfo%bsn_downstream(dsp+1:dsp+nthis) = thisinfo%bsn_downstream
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'basin_downstream', outinfo%bsn_downstream, &
            'catchment', compress = 1)

         ! ----- output : number of neighbours of a basin -----
         allocate(outinfo%bsn_num_nbr (outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) outinfo%bsn_num_nbr(dsp+1:dsp+nthis) = thisinfo%bsn_num_nbr
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'num_neighbour', outinfo%bsn_num_nbr, &
            'catchment', compress = 1)

         ! ----- output : neighbour index of a basin -----
         allocate(outinfo%bsn_idx_nbr (nnbmax, outinfo%ntotalcat))
         outinfo%bsn_idx_nbr(:,:) = -1
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) ndim1 = size(thisinfo%bsn_idx_nbr,1)
            IF (nthis > 0) outinfo%bsn_idx_nbr(1:ndim1,dsp+1:dsp+nthis) = thisinfo%bsn_idx_nbr
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'idx_neighbour', outinfo%bsn_idx_nbr, &
            'neighbour', 'catchment', compress = 1)

         ! ----- output : length of border between a basin and its neighbours -----
         allocate(outinfo%bsn_len_bdr (nnbmax, outinfo%ntotalcat))
         outinfo%bsn_len_bdr(:,:) = 0
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) ndim1 = size(thisinfo%bsn_len_bdr,1)
            IF (nthis > 0) outinfo%bsn_len_bdr(1:ndim1,dsp+1:dsp+nthis) = thisinfo%bsn_len_bdr
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'len_border', outinfo%bsn_len_bdr, &
            'neighbour', 'catchment', compress = 1)

         ! ----- output : basin elevation -----
         allocate(outinfo%bsn_elva (outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) outinfo%bsn_elva(dsp+1:dsp+nthis) = thisinfo%bsn_elva
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'basin_elevation', outinfo%bsn_elva, &
            'catchment', compress = 1)

         ! ----- output : lake id -----
         allocate(outinfo%lake_id (outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            IF (nthis > 0) outinfo%lake_id(dsp+1:dsp+nthis) = thisinfo%lake_id
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'lake_id', outinfo%lake_id, 'catchment', compress = 1)

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
      CALL mpi_finalize(p_err)

   END SUBROUTINE output_result

END MODULE output_mod
