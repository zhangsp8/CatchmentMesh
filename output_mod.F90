MODULE output_mod

CONTAINS
   
   SUBROUTINE output_result (maxhunum, maxnnb)

      USE task_mod
      USE hydro_data_mod
      USE ncio_serial

      IMPLICIT NONE

      integer (kind=4), intent(in) :: maxhunum
      integer (kind=4), intent(in) :: maxnnb

      integer (kind=4) :: np, mp
      real    (kind=8), allocatable :: longitude(:)
      real    (kind=8), allocatable :: latitude (:)
      integer (kind=4), allocatable :: catch (:,:)    
      real    (kind=4), allocatable :: hnd   (:,:)
      real    (kind=4), allocatable :: elv   (:,:)
      integer (kind=4), allocatable :: hunit (:,:)    

      type(info_typ) :: outinfo

      character(len=256) :: filename, mkdirname
      integer :: nthis, dsp, iunit, iblk, jblk, iloc, ic, ndim1
      logical :: found
      real (kind = 4), parameter :: spval = -1.e36

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) write(*,'(/A)') 'Step 5: Output results ...'

         
      IF (p_is_master) THEN
     
         IF (trim(storage_type) == 'one') THEN

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

            filename = trim(output_dir) // '/' // trim(casename) // '.nc'

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

            CALL ncio_define_dimension (filename, 'latitude', np)
            CALL ncio_define_dimension (filename, 'longitude', mp)
            CALL ncio_write_serial (filename, 'latitude', latitude,  'latitude')
            CALL ncio_write_serial (filename, 'longitude', longitude, 'longitude')
            
            CALL ncio_put_attr_str (filename, 'latitude', 'long_name', 'latitude')
            CALL ncio_put_attr_str (filename, 'latitude', 'units', 'degrees_north')
            CALL ncio_put_attr_str (filename, 'longitude', 'long_name', 'longitude')
            CALL ncio_put_attr_str (filename, 'longitude', 'units', 'degrees_east')


            CALL ncio_write_serial (filename, 'icatchment2d', catch, 'latitude', 'longitude', compress = 1)
            CALL ncio_write_serial (filename, 'hand',   hnd  , 'latitude', 'longitude', compress = 1)
            CALL ncio_write_serial (filename, 'elva',   elv  , 'latitude', 'longitude', compress = 1)
            CALL ncio_write_serial (filename, 'ihydrounit2d', hunit, 'latitude', 'longitude', compress = 1)

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

         filename = trim(output_dir) // '/' // trim(casename) // '.nc'
         IF (trim(storage_type) == 'block') THEN
            CALL ncio_create_file (filename)
         ENDIF

         CALL ncio_define_dimension (filename, 'catchment', outinfo%ntotalcat)
         CALL ncio_define_dimension (filename, 'hydrounit', maxhunum)
         CALL ncio_define_dimension (filename, 'neighbour', maxnnb)
         CALL ncio_define_dimension (filename, 'ncds', 2)


         ! ----- output : hydro unit index -----
         allocate(outinfo%hru_indx (maxhunum, outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            outinfo%hru_indx(:,dsp+1:dsp+nthis) = thisinfo%hru_indx
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_index', outinfo%hru_indx, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit area -----
         allocate(outinfo%hru_area (maxhunum, outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            outinfo%hru_area(:,dsp+1:dsp+nthis) = thisinfo%hru_area
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_area', outinfo%hru_area, &
            'hydrounit', 'catchment', compress = 1)
         
         ! ----- output : hydro unit hand -----
         allocate(outinfo%hru_hand (maxhunum, outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            outinfo%hru_hand(:,dsp+1:dsp+nthis) = thisinfo%hru_hand
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_hand', outinfo%hru_hand, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit elevation -----
         allocate(outinfo%hru_elva (maxhunum, outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            outinfo%hru_elva(:,dsp+1:dsp+nthis) = thisinfo%hru_elva
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_elva', outinfo%hru_elva, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit downstream -----
         allocate(outinfo%hru_next (maxhunum, outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            outinfo%hru_next(:,dsp+1:dsp+nthis) = thisinfo%hru_next
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_downstream', outinfo%hru_next, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit path length -----
         allocate(outinfo%hru_plen (maxhunum, outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            outinfo%hru_plen(:,dsp+1:dsp+nthis) = thisinfo%hru_plen
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_pathlen', outinfo%hru_plen, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : hydro unit drianage contour length -----
         allocate(outinfo%hru_lfac (maxhunum, outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            outinfo%hru_lfac(:,dsp+1:dsp+nthis) = thisinfo%hru_lfac
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'hydrounit_facelen', outinfo%hru_lfac, &
            'hydrounit', 'catchment', compress = 1)

         ! ----- output : river length -----
         allocate(outinfo%riv_len (outinfo%ntotalcat));   outinfo%riv_len(:) = spval
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%nrivseg
            outinfo%riv_len(dsp+1:dsp+nthis) = thisinfo%riv_len
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
            outinfo%riv_elv(dsp+1:dsp+nthis) = thisinfo%riv_elv
            dsp = dsp + thisinfo%ntotalcat
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'river_elevation', outinfo%riv_elv, &
            'catchment', compress = 1)
         CALL ncio_put_attr_real4 (filename, 'river_elevation', 'missing_value', spval)

         ! ----- output : number of hydro units in a basin -----
         allocate(outinfo%bsn_num_hru (outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            outinfo%bsn_num_hru(dsp+1:dsp+nthis) = thisinfo%bsn_num_hru
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'basin_numhru', outinfo%bsn_num_hru, &
            'catchment', compress = 1)

         ! ----- output : downstream basin -----
         allocate(outinfo%bsn_downstream (outinfo%ntotalcat))
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            outinfo%bsn_downstream(dsp+1:dsp+nthis) = thisinfo%bsn_downstream
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
            outinfo%bsn_num_nbr(dsp+1:dsp+nthis) = thisinfo%bsn_num_nbr
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'num_neighbour', outinfo%bsn_num_nbr, &
            'catchment', compress = 1)

         ! ----- output : neighbour index of a basin -----
         allocate(outinfo%bsn_idx_nbr (maxnnb, outinfo%ntotalcat))
         outinfo%bsn_idx_nbr(:,:) = -1
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            ndim1 = size(thisinfo%bsn_idx_nbr,1)
            outinfo%bsn_idx_nbr(1:ndim1,dsp+1:dsp+nthis) = thisinfo%bsn_idx_nbr
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO

         CALL ncio_write_serial (filename, 'idx_neighbour', outinfo%bsn_idx_nbr, &
            'neighbour', 'catchment', compress = 1)

         ! ----- output : length of border between a basin and its neighbours -----
         allocate(outinfo%bsn_len_bdr (maxnnb, outinfo%ntotalcat))
         outinfo%bsn_len_bdr(:,:) = 0
         thisinfo => allinfo;  dsp = 0
         DO WHILE (associated(thisinfo))
            nthis = thisinfo%ntotalcat
            ndim1 = size(thisinfo%bsn_len_bdr,1)
            outinfo%bsn_len_bdr(1:ndim1,dsp+1:dsp+nthis) = thisinfo%bsn_len_bdr
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
            outinfo%bsn_elva(dsp+1:dsp+nthis) = thisinfo%bsn_elva
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
            outinfo%lake_id(dsp+1:dsp+nthis) = thisinfo%lake_id
            dsp = dsp + nthis
            thisinfo => thisinfo%next
         ENDDO
         
         CALL ncio_write_serial (filename, 'lake_id', outinfo%lake_id, 'catchment', compress = 1)

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
      CALL mpi_finalize(p_err)

   END SUBROUTINE output_result 

END MODULE output_mod
