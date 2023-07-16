MODULE output_mod

CONTAINS
   
   SUBROUTINE output_result (output_dir, casename, storage_type, maxhunum, maxnnb)

      use task_mod
      use hydro_data_mod
      USE ncio_serial

      implicit none

      character(len=*), intent(in) :: output_dir, casename, storage_type
      integer (kind=4), intent(in) :: maxhunum
      integer (kind=4), intent(in) :: maxnnb

      integer (kind=4) :: np, mp
      REAL    (kind=8), allocatable :: longitude(:)
      REAL    (kind=8), allocatable :: latitude (:)
      integer (kind=4), allocatable :: catch (:,:)    
      real    (kind=4), allocatable :: hnd   (:,:)
      real    (kind=4), allocatable :: elv   (:,:)
      integer (kind=4), allocatable :: hunit (:,:)    

      character(len=256) :: filename, mkdirname
      integer :: iunit, iblk, jblk, iloc, nlat, nlon
      LOGICAL :: found

      if (p_is_master) write(*,'(/A)') 'Step 5: Output results ...'

      ! if (p_is_master) then
      !    filename = trim(output_dir) // '/' // trim(casename) // '_riv.bin'
      !    iunit = 100
      !    open  (iunit, file = trim(filename), form = 'unformatted', &
      !       status = 'replace', access = 'stream', action = 'write')
      !    write (iunit) riv_info_pix (:,1:nrivpix)
      !    close (iunit)
      ! end if

      if (trim(storage_type) == 'one') then

         ! shrink regions
         DO WHILE (.true.)
            iblk = (inorth-1)/nbox + 1
            found = .false.
            DO jblk = 1, mblock
               if (bkid(iblk,jblk) == p_iam_glb) then
                  iloc = inorth - blks(iblk,jblk)%idsp
                  IF (any(blks(iblk,jblk)%icat(iloc,:) > 0)) THEN
                     found = .true.
                     exit
                  ENDIF
               ENDIF
            ENDDO

            CALL mpi_allreduce (MPI_IN_PLACE, found, 1, MPI_LOGICAL, MPI_LOR, p_comm_glb, p_err)

            IF (found) THEN
               exit
            ELSE
               inorth = inorth + 1
            ENDIF
         ENDDO 

         DO WHILE (.true.)
            iblk = (isouth-1)/nbox + 1
            found = .false.
            DO jblk = 1, mblock
               if (bkid(iblk,jblk) == p_iam_glb) then
                  iloc = isouth - blks(iblk,jblk)%idsp
                  IF (any(blks(iblk,jblk)%icat(iloc,:) > 0)) THEN
                     found = .true.
                     exit
                  ENDIF
               ENDIF
            ENDDO

            CALL mpi_allreduce (MPI_IN_PLACE, found, 1, MPI_LOGICAL, MPI_LOR, p_comm_glb, p_err)

            IF (found) THEN
               exit
            ELSE
               isouth = isouth - 1
            ENDIF
         ENDDO 

         IF (jwest /= jeast) then
            DO WHILE (.true.)
               jblk = (jwest-1)/mbox + 1
               found = .false.
               DO iblk = 1, nblock
                  if (bkid(iblk,jblk) == p_iam_glb) then
                     iloc = jwest - blks(iblk,jblk)%jdsp
                     IF (any(blks(iblk,jblk)%icat(:,iloc) > 0)) THEN
                        found = .true.
                        exit
                     ENDIF
                  ENDIF
               ENDDO

               CALL mpi_allreduce (MPI_IN_PLACE, found, 1, MPI_LOGICAL, MPI_LOR, p_comm_glb, p_err)

               IF (found) THEN
                  exit
               ELSE
                  jwest = jwest + 1
                  IF (jwest == mglb+1) jwest = 1
               ENDIF
            ENDDO 

            DO WHILE (.true.)
               jblk = (jeast-1)/mbox + 1
               found = .false.
               DO iblk = 1, nblock
                  if (bkid(iblk,jblk) == p_iam_glb) then
                     iloc = jeast - blks(iblk,jblk)%jdsp
                     IF (any(blks(iblk,jblk)%icat(:,iloc) > 0)) THEN
                        found = .true.
                        exit
                     ENDIF
                  ENDIF
               ENDDO

               CALL mpi_allreduce (MPI_IN_PLACE, found, 1, MPI_LOGICAL, MPI_LOR, p_comm_glb, p_err)

               IF (found) THEN
                  exit
               ELSE
                  jeast = jeast - 1
                  IF (jeast == 0) jeast = mglb
               ENDIF
            ENDDO 
         ENDIF

         if (p_is_master) then

            np = isouth - inorth + 1
            if (jeast == jwest) then
               mp = mglb
            else
               mp = jeast  - jwest  + 1
               if (mp < 0) mp = mp + mglb 
            end if

            allocate (latitude (np))
            allocate (longitude(mp))
            allocate (catch (np,mp))
            allocate (hnd   (np,mp))
            allocate (elv   (np,mp))
            allocate (hunit (np,mp))

            if (jwest /= jeast) then
               call aggregate_data (inorth, isouth, jwest, jeast, & 
                  np, mp, longitude, latitude, &
                  icat = catch, hnd = hnd, elv = elv, hunit = hunit)
            else
               call aggregate_data (inorth, isouth, 1, mglb, & 
                  np, mp, longitude, latitude, &
                  icat = catch, hnd = hnd, elv = elv, hunit = hunit)
            end if

            filename = trim(output_dir) // '/' // trim(casename) // '.nc'

            call ncio_create_file (filename)
            call ncio_write_serial (filename, 'inorth', inorth)
            call ncio_write_serial (filename, 'isouth', isouth)
            call ncio_write_serial (filename, 'jwest', jwest)
            call ncio_write_serial (filename, 'jeast', jeast)

            nlat = isouth - inorth + 1
            nlon = jeast - jwest + 1
            IF (nlon < 0) nlon = nlon + mglb
            CALL ncio_define_dimension (filename, 'latitude',  nlat)
            CALL ncio_define_dimension (filename, 'longitude', nlon)
            call ncio_write_serial (filename, 'latitude',  latitude,  'latitude')
            call ncio_write_serial (filename, 'longitude', longitude, 'longitude')

            call ncio_write_serial (filename, 'icatchment2d', catch, 'latitude', 'longitude', compress = 1)
            call ncio_write_serial (filename, 'hand',   hnd  , 'latitude', 'longitude', compress = 1)
            call ncio_write_serial (filename, 'elva',   elv  , 'latitude', 'longitude', compress = 1)
            call ncio_write_serial (filename, 'ihydrounit2d', hunit, 'latitude', 'longitude', compress = 1)

            deallocate (latitude )
            deallocate (longitude)
            deallocate (catch)
            deallocate (hnd  )
            deallocate (elv  )
            deallocate (hunit)

            call excute_data_task (t_exit)

         elseif (p_is_data) then
            call data_daemon ()
         end if

      elseif (trim(storage_type) == 'block') then

         if (p_is_master) then

            mkdirname = 'mkdir -p ' // trim(output_dir) // '/' // trim(casename)
            call system(trim(mkdirname))

            call mpi_barrier (p_comm_data, p_err)

         elseif (p_is_data) then

            call mpi_barrier (p_comm_data, p_err)

            do jblk = 1, mblock
               do iblk = 1, nblock
                  if (bkid(iblk,jblk) == p_iam_glb) then

                     call get_filename (trim(output_dir) // '/' // trim(casename), &
                        iblk, jblk, filename)

                     call ncio_create_file (filename)
                     CALL ncio_define_dimension (filename, 'latitude',  nbox)
                     CALL ncio_define_dimension (filename, 'longitude', mbox)
                     call ncio_write_serial (filename, 'latitude',  blks(iblk,jblk)%lat, 'latitude' ) 
                     call ncio_write_serial (filename, 'longitude', blks(iblk,jblk)%lon, 'longitude')

                     call ncio_write_serial (filename, 'icatchment2d', blks(iblk,jblk)%icat, &
                        'latitude', 'longitude', compress = 1)
                     call ncio_write_serial (filename, 'hand',   blks(iblk,jblk)%hnd,  'latitude', 'longitude', compress = 1)
                     call ncio_write_serial (filename, 'elva',   blks(iblk,jblk)%elv,  'latitude', 'longitude', compress = 1)
                     call ncio_write_serial (filename, 'ihydrounit2d', blks(iblk,jblk)%hunit, &
                        'latitude', 'longitude', compress = 1)

                  end if
               end do
            end do

         end if

      end if

      IF (p_is_master) THEN

         filename = trim(output_dir) // '/' // trim(casename) // '.nc'
         if (trim(storage_type) == 'block') then
            CALL ncio_create_file (filename)
         ENDIF

         CALL ncio_define_dimension (filename, 'catchment', ntotalcat)
         CALL ncio_define_dimension (filename, 'hydrounit', maxhunum)
         CALL ncio_define_dimension (filename, 'neighbour', maxnnb)
         
         call ncio_write_serial (filename, 'lake_id', lake_info_id, 'catchment', compress = 1)

         call ncio_write_serial (filename, 'hydrounit_index', hru_info_indx, &
            'hydrounit', 'catchment', compress = 1)
         call ncio_write_serial (filename, 'hydrounit_area', hru_info_area, &
            'hydrounit', 'catchment', compress = 1)
         call ncio_write_serial (filename, 'hydrounit_downstream', hru_info_next, &
            'hydrounit', 'catchment', compress = 1)
         call ncio_write_serial (filename, 'hydrounit_hand', hru_info_hand, &
            'hydrounit', 'catchment', compress = 1)
         call ncio_write_serial (filename, 'hydrounit_elva', hru_info_elva, &
            'hydrounit', 'catchment', compress = 1)
         call ncio_write_serial (filename, 'hydrounit_pathlen', hru_info_plen, &
            'hydrounit', 'catchment', compress = 1)
         call ncio_write_serial (filename, 'hydrounit_facelen', hru_info_lfac, &
            'hydrounit', 'catchment', compress = 1)

         call ncio_write_serial (filename, 'river_length',    riv_info_len, 'catchment', compress = 1)
         call ncio_write_serial (filename, 'river_elevation', riv_info_elv, 'catchment', compress = 1)

         call ncio_write_serial (filename, 'basin_numhru', bsn_info_num_hru, 'catchment', compress = 1)
         call ncio_write_serial (filename, 'basin_downstream', bsn_info_downstream, 'catchment', compress = 1)
         call ncio_write_serial (filename, 'basin_num_neighbour', bsn_info_num_nbr, 'catchment', compress = 1)
         call ncio_write_serial (filename, 'basin_idx_neighbour', bsn_info_idx_nbr, &
            'neighbour', 'catchment', compress = 1)
         call ncio_write_serial (filename, 'basin_len_border', bsn_info_len_bdr, &
            'neighbour', 'catchment', compress = 1)

      ENDIF

   END SUBROUTINE output_result 

END MODULE output_mod
