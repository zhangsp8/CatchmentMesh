MODULE hydro_data_mod

   USE ncio_serial
   IMPLICIT NONE

   integer (kind=4), parameter :: nbox = 6000
   integer (kind=4), parameter :: mbox = 6000
   integer (kind=4), parameter :: nglb = 216000
   integer (kind=4), parameter :: mglb = 432000
   integer (kind=4), parameter :: nblock = 36
   integer (kind=4), parameter :: mblock = 72

   character(len=256) :: def_hydro_dir,  def_lake_dir
   character(len=256) :: def_output_dir, def_casename, def_storage_type

   type block_typ

      integer :: idsp, jdsp
      logical :: ready, is_ocean

      real (kind=8), allocatable :: lat(:)
      real (kind=8), allocatable :: lon(:)

      ! Flow direction is prepared in 1-byte SIGNED integer (int8) and is defined as:
      !    1: east, 2: southeast, 4: south, 8: southwest, 16: west, 32: northwest, 64: north. 128: northeast
      !    0: river mouth, -1: inland depression, -9: undefined (ocean)
      ! IF a flow direction file is opened as UNSIGNED integer, undefined=247 and inland depression=255
      integer(kind=1), allocatable :: dir  (:,:)  ! flow direction
      real   (kind=4), allocatable :: upa  (:,:)  ! upstream drainage area (km^2)
      real   (kind=4), allocatable :: elv  (:,:)  ! elevation (m)
      real   (kind=4), allocatable :: wth  (:,:)  ! river width (m)
      integer(kind=4), allocatable :: lake (:,:)  ! lake id

      ! normal catchment with upstream area larger than maximum catchment size (>0);
      ! flow to ocean(-9); flow to inland depression(-1);
      ! flow to small river mouth (-2); flow outside the region(-3)
      integer(kind=4), allocatable :: icat (:,:)  ! catchment number
      real   (kind=4), allocatable :: hnd  (:,:)  ! height above nearest drainage (m)
      integer(kind=4), allocatable :: hunit(:,:)  ! hillslope units

   END type block_typ


   type (block_typ) :: blks (nblock,mblock)

   integer (kind=4) :: inorth, isouth, jwest, jeast
   logical :: def_include_all_upstream = .false.

   real (kind=4) :: dlat(nglb)
   real (kind=4) :: dlon(nglb)

   type info_typ

      integer :: ithisblk, jthisblk

      ! riv_info_pix(1:2,:) : river pixel coordinate;
      ! riv_info_pix(3,:)   : catchment id this pixel located in
      integer (kind=4) :: nrivpix, nrivseg
      integer (kind=4), allocatable :: riv_pix (:,:)
      real    (kind=4), allocatable :: riv_len (:)
      real    (kind=4), allocatable :: riv_elv (:)
      real    (kind=4), allocatable :: riv_dep (:)

      ! hydrounit information
      integer (kind=4), allocatable :: bsn_num_hru(:)
      integer (kind=4), allocatable :: hru_indx (:,:)
      real    (kind=4), allocatable :: hru_area (:,:)
      real    (kind=4), allocatable :: hru_hand (:,:)
      real    (kind=4), allocatable :: hru_elva (:,:)
      integer (kind=4), allocatable :: hru_next (:,:)
      real    (kind=4), allocatable :: hru_plen (:,:)
      real    (kind=4), allocatable :: hru_lfac (:,:)

      ! basin information
      integer (kind=4) :: icatdsp
      integer (kind=4) :: ntotalcat
      integer (kind=4), allocatable :: bsn_index      (:)
      integer (kind=4), allocatable :: bsn_nswe     (:,:)
      integer (kind=4), allocatable :: bsn_downstream (:)
      integer (kind=4), allocatable :: bsn_num_nbr    (:)
      integer (kind=4), allocatable :: bsn_idx_nbr  (:,:)
      real    (kind=4), allocatable :: bsn_len_bdr  (:,:)
      real    (kind=4), allocatable :: bsn_elva       (:)

      ! lake information
      integer (kind=4), allocatable :: lake_id (:)

      type(info_typ), pointer :: next

   END type info_typ

   type(info_typ), target  :: allinfo
   type(info_typ), pointer :: thisinfo

CONTAINS

   !----------------------------------------
   SUBROUTINE get_filename (dirname, iblk, jblk, filename)

      IMPLICIT NONE

      integer (kind=4), intent(in)  :: iblk, jblk
      character(len=*), intent(in)  :: dirname
      character(len=*), intent(out) :: filename
      character(len=3)   :: pre1
      character(len=4)   :: pre2

      IF (iblk <= 18) THEN
         write (pre1,'(A1,I2.2)') 'n', (18-iblk)*5
      ELSE
         write (pre1,'(A1,I2.2)') 's', (iblk-18)*5
      ENDIF

      IF (jblk <= 36) THEN
         write (pre2,'(A1,I3.3)') 'w', (37-jblk)*5
      ELSE
         write (pre2,'(A1,I3.3)') 'e', (jblk-37)*5
      ENDIF

      filename = trim(dirname) // '/' // trim(pre1) // trim(pre2) // '.nc'


   END SUBROUTINE get_filename

   !----------------------------------------
   SUBROUTINE get_region (west, east, north, south)

      USE task_mod
      IMPLICIT NONE

      real (kind=4), intent(in) :: west,  east    ! -180 to 180
      real (kind=4), intent(in) :: north, south   ! -90 to 90

      integer (kind=4) :: iblk, jblk

      IF (p_is_master) THEN

         IF (west == east) THEN
            jwest = 1
            jeast = mglb
         ELSE
            jwest = floor((west + 180.0004167) * 1200) + 1
            jeast = ceiling((east + 180.0004167) * 1200)

            jwest = max(1,min(jwest,mglb))
            jeast = max(1,min(jeast,mglb))
         ENDIF

         inorth = floor((89.9995833 - north) * 1200) + 1
         isouth = ceiling((89.9995833 - south) * 1200)

         inorth = max(1,min(inorth,nglb))
         isouth = max(1,min(isouth,nglb))

         write(*,'(A,4F8.2)') 'Step 0 : Region (n,s,w,e): ', north, south, west, east
         write(*,'(A,4I8)')   'Step 0 : Iregion(n,s,w,e): ', inorth, isouth, jwest, jeast

         DO jblk = 1, mblock
            DO iblk = 1, nblock
               blks(iblk,jblk)%idsp = (iblk-1)*nbox
               blks(iblk,jblk)%jdsp = (jblk-1)*mbox
               blks(iblk,jblk)%ready = .false.
            ENDDO
         ENDDO

      ENDIF

   END SUBROUTINE get_region

   !----------------------------------------
   logical FUNCTION on_region_boundary (i, j)

      IMPLICIT NONE
      integer :: i, j

      IF (within_region(i,j,.true.)) THEN
         on_region_boundary = &
            ((i == inorth) .and. (inorth /= 1))   &
            .or. ((i == isouth) .and. (isouth /= 180000))   &
            .or. (((j == jwest) .or. (j == jeast )) .and. ((jwest /= 1) .or. (jeast /= mglb)))
      ELSE
         on_region_boundary = .false.
      ENDIF

   END FUNCTION

   !----------------------------------------
   SUBROUTINE readin_blocks (iblk, jblk)

      USE task_mod
      IMPLICIT NONE

      integer (kind=4) :: iproc, iblk, jblk, i, j
      integer (kind=4) :: iwin, jwin

      character(len=256) :: filename, lakefile
      logical :: fexists, load_from_file


      CALL get_filename (def_hydro_dir, iblk, jblk, filename)

      inquire (file=trim(filename), exist=fexists)
      IF (fexists) THEN

         CALL ncio_read_serial (filename, 'latitude',  blks(iblk,jblk)%lat)
         CALL ncio_read_serial (filename, 'longitude', blks(iblk,jblk)%lon)

         CALL ncio_read_serial (filename, 'dir', blks(iblk,jblk)%dir)
         CALL ncio_read_serial (filename, 'upa', blks(iblk,jblk)%upa)
         CALL ncio_read_serial (filename, 'elv', blks(iblk,jblk)%elv)
         CALL ncio_read_serial (filename, 'wth', blks(iblk,jblk)%wth)

         CALL get_filename (def_lake_dir, iblk, jblk, lakefile)
         CALL ncio_read_serial (lakefile, 'lake', blks(iblk,jblk)%lake)

      ELSE

         allocate (blks(iblk,jblk)%lat(mbox));  blks(iblk,jblk)%lat = -1.e36
         allocate (blks(iblk,jblk)%lon(mbox));  blks(iblk,jblk)%lon = -1.e36

         allocate (blks(iblk,jblk)%dir  (nbox,mbox)); blks(iblk,jblk)%dir  = -9
         allocate (blks(iblk,jblk)%upa  (nbox,mbox)); blks(iblk,jblk)%upa  = 0.
         allocate (blks(iblk,jblk)%elv  (nbox,mbox)); blks(iblk,jblk)%elv  = 0.
         allocate (blks(iblk,jblk)%wth  (nbox,mbox)); blks(iblk,jblk)%wth  = 0.
         allocate (blks(iblk,jblk)%lake (nbox,mbox)); blks(iblk,jblk)%lake = 0

      ENDIF

      blks(iblk,jblk)%is_ocean = (.not. fexists)

      IF (trim(def_storage_type) == 'block') THEN
         CALL get_filename (trim(def_output_dir) // '/' // trim(def_casename), iblk, jblk, filename)
         inquire (file=trim(filename), exist=load_from_file)
      ELSE
         load_from_file = .false.
      ENDIF

      IF (load_from_file) THEN

         CALL ncio_read_serial (filename, 'icatchment2d', blks(iblk,jblk)%icat )
         CALL ncio_read_serial (filename, 'hand',         blks(iblk,jblk)%hnd  )
         CALL ncio_read_serial (filename, 'ihydrounit2d', blks(iblk,jblk)%hunit)

      ELSE

         allocate (blks(iblk,jblk)%hnd  (nbox,mbox))
         allocate (blks(iblk,jblk)%icat (nbox,mbox))
         allocate (blks(iblk,jblk)%hunit(nbox,mbox))

         IF (def_include_all_upstream) THEN
            blks(iblk,jblk)%icat  = 0
         ELSE
            blks(iblk,jblk)%icat  = -3
         ENDIF

         DO i = 1, nbox
            IF ((blks(iblk,jblk)%idsp + i >= inorth) &
               .and. (blks(iblk,jblk)%idsp + i <= isouth)) THEN

               DO j = 1, mbox
                  IF ((jwest > jeast) .and. (blks(iblk,jblk)%jdsp + j < jwest) &
                     .and. (blks(iblk,jblk)%jdsp + j > jeast)) THEN
                     CYCLE
                  ENDIF

                  IF ((jwest < jeast) .and. ((blks(iblk,jblk)%jdsp + j < jwest) &
                     .or. (blks(iblk,jblk)%jdsp + j > jeast))) THEN
                     CYCLE
                  ENDIF

                  IF (fexists) THEN
                     IF (blks(iblk,jblk)%dir(i,j) == -9) THEN
                        ! ocean
                        blks(iblk,jblk)%icat(i,j) = -9
                     ELSE
                        blks(iblk,jblk)%icat(i,j) = 0
                     ENDIF
                  ELSE
                     ! ocean
                     blks(iblk,jblk)%icat(i,j) = -9
                  ENDIF

               ENDDO
            ENDIF
         ENDDO

         blks(iblk,jblk)%hnd   = 0
         blks(iblk,jblk)%hunit = -1

      ENDIF

      blks(iblk,jblk)%ready = .true.

   END SUBROUTINE readin_blocks

   !----------------------------------------
   SUBROUTINE flush_blocks (output)

      IMPLICIT NONE

      logical, intent(in) :: output

      integer (kind=4) :: iblk, jblk
      logical :: fexists
      character(len=256) :: mkdirname, filename

      mkdirname = 'mkdir -p ' // trim(def_output_dir) // '/' // trim(def_casename)
      CALL system(trim(mkdirname))

      DO jblk = 1, mblock
         DO iblk = 1, nblock
            IF (blks(iblk,jblk)%ready) THEN

               IF (output) THEN

                  CALL get_filename (trim(def_output_dir) // '/' // trim(def_casename), &
                     iblk, jblk, filename)

                  write(*,*) 'Write cache data to file ', trim(filename)

                  inquire (file=trim(filename), exist=fexists)

                  IF (.not. fexists) THEN
                     CALL ncio_create_file (filename)
                     CALL ncio_define_dimension (filename, 'latitude',  nbox)
                     CALL ncio_define_dimension (filename, 'longitude', mbox)
                     CALL ncio_write_serial (filename, 'latitude',  blks(iblk,jblk)%lat, &
                        'latitude' )
                     CALL ncio_write_serial (filename, 'longitude', blks(iblk,jblk)%lon, &
                        'longitude')
                     CALL ncio_write_serial (filename, 'elva', blks(iblk,jblk)%elv,  &
                        'latitude', 'longitude', compress = 1)
                  ENDIF

                  CALL ncio_write_serial (filename, 'icatchment2d', blks(iblk,jblk)%icat, &
                     'latitude', 'longitude', compress = 1)
                  CALL ncio_write_serial (filename, 'hand', blks(iblk,jblk)%hnd, &
                     'latitude', 'longitude', compress = 1)
                  CALL ncio_write_serial (filename, 'ihydrounit2d', blks(iblk,jblk)%hunit, &
                     'latitude', 'longitude', compress = 1)

               ENDIF

               deallocate (blks(iblk,jblk)%lat)
               deallocate (blks(iblk,jblk)%lon)

               deallocate (blks(iblk,jblk)%dir )
               deallocate (blks(iblk,jblk)%upa )
               deallocate (blks(iblk,jblk)%elv )
               deallocate (blks(iblk,jblk)%wth )

               deallocate (blks(iblk,jblk)%hnd  )
               deallocate (blks(iblk,jblk)%icat )
               deallocate (blks(iblk,jblk)%hunit)

               deallocate (blks(iblk,jblk)%lake )

               blks(iblk,jblk)%ready = .false.

            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE flush_blocks

   !------------------------------------------------
   SUBROUTINE init_area ()

      USE MathConstants, only : deg2rad
      IMPLICIT NONE

      real(kind=8), parameter :: re = 6.37122e3 ! kilometer
      real(kind=8) :: lat0, lat1
      integer (kind=4) :: ilat

      DO ilat = 1, nglb
         lat1 = 90.0 - 1 / 1200.0 / 2.0 - 180.0 / nglb *  (ilat-1)
         lat0 = lat1 - 1 / 1200.0
         lat1 = deg2rad * min( 90.0, lat1)
         lat0 = deg2rad * max(-90.0, lat0)

         dlat(ilat) = re * (lat1 - lat0)
         dlon(ilat) = re * (sin(lat1)-sin(lat0)) / (lat1-lat0) * deg2rad/1200.0
      ENDDO

   END SUBROUTINE init_area

   !-------------------------------------------------
   SUBROUTINE nextij (i, j, dir, inext, jnext, is_local)

      IMPLICIT NONE

      integer(kind=4), intent(in) :: i, j
      integer(kind=1), intent(in) :: dir

      integer(kind=4), intent(out)  :: inext, jnext

      logical, optional, intent(in) :: is_local

      select CASE (dir)
      CASE (1)
         inext = i
         jnext = j + 1
      CASE (2)
         inext = i + 1
         jnext = j + 1
      CASE (4)
         inext = i + 1
         jnext = j
      CASE (8)
         inext = i + 1
         jnext = j - 1
      CASE (16)
         inext = i
         jnext = j - 1
      CASE (32)
         inext = i - 1
         jnext = j - 1
      CASE (64)
         inext = i - 1
         jnext = j
      CASE (-128)
         inext = i - 1
         jnext = j + 1
      CASE default
         inext = i
         jnext = j
      END select

      IF (.not. present(is_local)) THEN
         IF (jnext == 0)   jnext = mglb
         IF (jnext > mglb) jnext = 1
      ELSEIF (.not. is_local) THEN
         IF (jnext == 0)   jnext = mglb
         IF (jnext > mglb) jnext = 1
      ENDIF

   END SUBROUTINE nextij

   !----------------------------------------------
   SUBROUTINE enlarge_nswe (imin, imax, jmin, jmax)

      IMPLICIT NONE

      integer, intent(inout) :: imin, imax, jmin, jmax

      IF (def_include_all_upstream) THEN
         imin = max(imin-1, 1)
         imax = min(imax+1, nglb)
      ELSE
         imin = max(imin-1, inorth)
         imax = min(imax+1, isouth)
      ENDIF

      IF ((jwest == 1) .and. (jeast == mglb)) THEN
         jmin = jmin-1
         jmax = jmax+1
      ELSE
         IF (def_include_all_upstream) THEN
            jmin = jmin-1
            jmax = jmax+1
         ELSE
            IF (jmin /= jwest) jmin = jmin-1
            IF (jmax /= jeast) jmax = jmax+1
         ENDIF
      ENDIF

      IF (jmin == 0)   jmin = mglb
      IF (jmax > mglb) jmax = 1

   END SUBROUTINE enlarge_nswe

   !----------------------------------------------
   logical FUNCTION is_feasible_step (i, j, i_dn, j_dn)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i, j, i_dn, j_dn

      IF ((i /= i_dn) .or. (j /= j_dn)) THEN
         is_feasible_step = .true.
      ELSE
         is_feasible_step = .false.
         RETURN
      ENDIF

      IF ((i_dn < inorth) .or. (i_dn > isouth)) THEN
         is_feasible_step = .false.
         RETURN
      ENDIF

      IF (jwest /= 0) THEN
         IF (jwest < jeast) THEN
            IF ((j_dn < jwest) .or. (j_dn > jeast)) THEN
               is_feasible_step = .false.
               RETURN
            ENDIF
         ELSE
            IF ((j_dn < jwest) .and. (j_dn > jeast)) THEN
               is_feasible_step = .false.
               RETURN
            ENDIF
         ENDIF
      ENDIF

   END FUNCTION is_feasible_step

   !----------------------------------------------
   logical FUNCTION within_region (i, j, strict)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i, j
      logical :: strict

      within_region = .true.

      IF ((.not. strict) .and. def_include_all_upstream) THEN
         RETURN
      ENDIF

      IF ((i < inorth) .or. (i > isouth)) THEN
         within_region = .false.
         RETURN
      ENDIF

      IF (jwest /= 0) THEN
         IF (jwest < jeast) THEN
            IF ((j < jwest) .or. (j > jeast)) THEN
               within_region = .false.
               RETURN
            ENDIF
         ELSE
            IF ((j < jwest) .and. (j > jeast)) THEN
               within_region = .false.
               RETURN
            ENDIF
         ENDIF
      ENDIF

   END FUNCTION within_region

   !---------------------------------------------
   SUBROUTINE block_iterator ( &
         imin, imax, jmin, jmax, iblk, jblk, &
         i0, i1, j0, j1, il0, il1, jl0, jl1, &
         iblk1, jblk1, end_of_data)

      IMPLICIT NONE

      integer, intent(in)  :: imin, imax, jmin, jmax
      integer, intent(in)  :: iblk, jblk
      integer, intent(out) :: i0, i1, j0, j1
      integer, intent(out) :: il0, il1, jl0, jl1
      integer, intent(out) :: iblk1, jblk1
      logical, intent(out) :: end_of_data

      i0 = max(imin-blks(iblk,jblk)%idsp, 1)
      i1 = min(imax-blks(iblk,jblk)%idsp, nbox)
      il0 = i0 + blks(iblk,jblk)%idsp - imin + 1
      il1 = i1 + blks(iblk,jblk)%idsp - imin + 1

      IF (jmax > blks(iblk,jblk)%jdsp) THEN
         j1 = min(jmax-blks(iblk,jblk)%jdsp, mbox)
      ELSE
         j1 = mbox
      ENDIF

      IF (jmin <= blks(iblk,jblk)%jdsp+mbox) THEN
         j0  = max(jmin-blks(iblk,jblk)%jdsp, 1)
         jl0 = j0 + blks(iblk,jblk)%jdsp - jmin + 1
         jl1 = j1 + blks(iblk,jblk)%jdsp - jmin + 1
      ELSE
         j0  = 1
         jl0 = j0 + blks(iblk,jblk)%jdsp + mglb - jmin + 1
         jl1 = j1 + blks(iblk,jblk)%jdsp + mglb - jmin + 1
      ENDIF

      iblk1 = iblk
      end_of_data = .false.
      IF (blks(iblk,jblk)%jdsp > jmax) THEN
         jblk1 = mod(jblk,mblock) + 1
      ELSEIF (blks(iblk,jblk)%jdsp+mbox < jmax) THEN
         jblk1 = jblk + 1
      ELSEIF (blks(iblk,jblk)%idsp+nbox < imax) THEN
         iblk1 = iblk1 + 1
         jblk1 = (jmin-1)/mbox + 1
      ELSE
         end_of_data = .true.
      ENDIF

   END SUBROUTINE block_iterator

   !----------------------------------------------
   integer(kind=1) FUNCTION get_dir (i, j)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iblk, jblk

      iblk = (i - 1)/nbox + 1
      jblk = (j - 1)/mbox + 1

      IF (.not. blks(iblk,jblk)%ready)  CALL readin_blocks(iblk, jblk)

      get_dir = blks(iblk,jblk)%dir (i-blks(iblk,jblk)%idsp, j-blks(iblk,jblk)%jdsp)

   END FUNCTION get_dir

   !----------------------------------------------
   real(kind=4) FUNCTION get_upa (i, j)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iblk, jblk

      iblk = (i - 1)/nbox + 1
      jblk = (j - 1)/mbox + 1

      IF (.not. blks(iblk,jblk)%ready)  CALL readin_blocks(iblk, jblk)

      get_upa = blks(iblk,jblk)%upa (i-blks(iblk,jblk)%idsp, j-blks(iblk,jblk)%jdsp)

   END FUNCTION get_upa

   !----------------------------------------------
   real(kind=4) FUNCTION get_elv (i, j)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iblk, jblk

      iblk = (i - 1)/nbox + 1
      jblk = (j - 1)/mbox + 1

      IF (.not. blks(iblk,jblk)%ready)  CALL readin_blocks(iblk, jblk)

      get_elv = blks(iblk,jblk)%elv (i-blks(iblk,jblk)%idsp, j-blks(iblk,jblk)%jdsp)

   END FUNCTION get_elv

   !----------------------------------------------
   integer (kind=4) FUNCTION get_icat (i, j)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iblk, jblk

      iblk = (i - 1)/nbox + 1
      jblk = (j - 1)/mbox + 1

      IF (.not. blks(iblk,jblk)%ready)  CALL readin_blocks(iblk, jblk)

      get_icat = blks(iblk,jblk)%icat (i-blks(iblk,jblk)%idsp, j-blks(iblk,jblk)%jdsp)

   END FUNCTION get_icat

   !----------------------------------------------
   real(kind=4) FUNCTION get_wth (i, j)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iblk, jblk

      iblk = (i - 1)/nbox + 1
      jblk = (j - 1)/mbox + 1

      IF (.not. blks(iblk,jblk)%ready)  CALL readin_blocks(iblk, jblk)

      get_wth = blks(iblk,jblk)%wth (i-blks(iblk,jblk)%idsp, j-blks(iblk,jblk)%jdsp)

   END FUNCTION get_wth

   !----------------------------------------------
   integer(kind=4) FUNCTION get_lake (i, j)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iblk, jblk

      iblk = (i - 1)/nbox + 1
      jblk = (j - 1)/mbox + 1

      IF (.not. blks(iblk,jblk)%ready)  CALL readin_blocks(iblk, jblk)

      get_lake = blks(iblk,jblk)%lake (i-blks(iblk,jblk)%idsp, j-blks(iblk,jblk)%jdsp)

   END FUNCTION get_lake

   !----------------------------------------------
   real(kind=4) FUNCTION get_area (i, j)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i, j

      get_area = dlon(i) * dlat(i)

   END FUNCTION get_area

   !----------------------------------------------
   real(kind=4) FUNCTION dist_between (i0, j0, i1, j1)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i0, j0, i1, j1

      IF (i0 == i1) THEN
         dist_between = dlon(i0)
      ELSEIF (j0 == j1) THEN
         dist_between = 0.5*(dlat(i0) + dlat(i1))
      ELSE
         dist_between = sqrt((0.5*(dlat(i0)+dlat(i1)))**2 + (0.5*(dlon(i0)+dlon(i1)))**2)
      ENDIF

   END FUNCTION dist_between

   !------------------------------------
   SUBROUTINE aggregate_data (imin, imax, jmin, jmax, &
         np, mp, longitude, latitude, icat, dir, hnd, elv, hunit)

      USE task_mod
      IMPLICIT NONE

      integer (kind=4), intent(in) :: imin, imax, jmin, jmax

      integer (kind=4), intent(in) :: np, mp

      real    (kind=8), intent(inout), optional :: longitude(mp)
      real    (kind=8), intent(inout), optional :: latitude (np)
      integer (kind=4), intent(inout), optional :: icat  (np,mp)
      integer (kind=1), intent(inout), optional :: dir   (np,mp)
      real    (kind=4), intent(inout), optional :: hnd   (np,mp)
      real    (kind=4), intent(inout), optional :: elv   (np,mp)
      integer (kind=4), intent(inout), optional :: hunit (np,mp)

      integer (kind=4) :: iblk, jblk, iblk1, jblk1
      integer (kind=4) :: i0, i1, j0, j1, il0, il1, jl0, jl1
      logical :: end_of_data

      IF (present(icat)) THEN
         icat(:,:) = -9
      ENDIF

      iblk = (imin-1)/nbox + 1
      jblk = (jmin-1)/mbox + 1
      DO WHILE (.true.)

         CALL block_iterator (imin, imax, jmin, jmax, iblk, jblk, &
            i0, i1, j0, j1, il0, il1, jl0, jl1, &
            iblk1, jblk1, end_of_data)

         IF (.not. blks(iblk,jblk)%ready)  CALL readin_blocks(iblk, jblk)

         IF (present(longitude)) then
            IF (.not. any(blks(iblk,jblk)%lon(j0:j1) == -1.e36)) THEN
               longitude(jl0:jl1) = blks(iblk,jblk)%lon(j0:j1)
            ENDIF
         ENDIF

         IF (present(latitude )) then
            IF (.not. any(blks(iblk,jblk)%lat(i0:i1) == -1.e36)) THEN
               latitude (il0:il1) = blks(iblk,jblk)%lat(i0:i1)
            ENDIF
         ENDIF

         IF (present(icat))  icat (il0:il1,jl0:jl1) = blks(iblk,jblk)%icat (i0:i1,j0:j1)
         IF (present(dir))   dir  (il0:il1,jl0:jl1) = blks(iblk,jblk)%dir  (i0:i1,j0:j1)
         IF (present(hnd))   hnd  (il0:il1,jl0:jl1) = blks(iblk,jblk)%hnd  (i0:i1,j0:j1)
         IF (present(elv))   elv  (il0:il1,jl0:jl1) = blks(iblk,jblk)%elv  (i0:i1,j0:j1)
         IF (present(hunit)) hunit(il0:il1,jl0:jl1) = blks(iblk,jblk)%hunit(i0:i1,j0:j1)

         IF (end_of_data) THEN
            EXIT
         ELSE
            iblk = iblk1; jblk = jblk1
         ENDIF
      ENDDO

   END SUBROUTINE aggregate_data

   !------------------------------------
   SUBROUTINE aggregate_catch (imin, imax, jmin, jmax, np, mp, icat)

      USE task_mod
      IMPLICIT NONE

      integer (kind=4), intent(in) :: imin, imax, jmin, jmax
      integer (kind=4), intent(in) :: np, mp
      integer (kind=4), intent(inout) :: icat(np,mp)

      integer (kind=4) :: iblk, jblk, iblk1, jblk1
      integer (kind=4) :: i0, i1, j0, j1, il0, il1, jl0, jl1
      logical :: end_of_data
      character(len=256) :: filename

      icat(:,:) = -9

      iblk = (imin-1)/nbox + 1
      jblk = (jmin-1)/mbox + 1
      DO WHILE (.true.)

         CALL block_iterator (imin, imax, jmin, jmax, iblk, jblk, &
            i0, i1, j0, j1, il0, il1, jl0, jl1, &
            iblk1, jblk1, end_of_data)

         IF (.not. allocated(blks(iblk,jblk)%icat)) THEN
            CALL get_filename (trim(def_output_dir) // '/' // trim(def_casename), iblk, jblk, filename)
            CALL ncio_read_serial (filename, 'icatchment2d', blks(iblk,jblk)%icat )
         ENDIF

         icat(il0:il1,jl0:jl1) = blks(iblk,jblk)%icat(i0:i1,j0:j1)

         IF (end_of_data) THEN
            EXIT
         ELSE
            iblk = iblk1; jblk = jblk1
         ENDIF
      ENDDO

   END SUBROUTINE aggregate_catch

   !-------------------------------------------------------
   SUBROUTINE append_river (river, nriv, ij0, icatch)

      IMPLICIT NONE

      integer (kind=4), intent(inout), allocatable :: river(:,:)
      integer (kind=4), intent(inout) :: nriv
      integer (kind=4), intent(in) :: ij0(2), icatch

      integer (kind=4) :: nsize
      integer (kind=4), allocatable :: temp(:,:)

      nriv = nriv + 1

      nsize = size(river,2)
      IF (nriv > nsize) THEN
         allocate (temp(3,nsize))
         temp = river

         deallocate (river)
         allocate   (river(3,nsize+100000))
         river(:,1:nsize) = temp

         deallocate (temp)
      ENDIF

      river (1:2,nriv) = ij0
      river   (3,nriv) = icatch

   END SUBROUTINE append_river

   !--------------------------------------------------------------
   SUBROUTINE free_memory ()

      USE task_mod
      IMPLICIT NONE

      integer (kind=4) :: iblk, jblk


      DO iblk = 1, nblock
         DO jblk = 1, mblock
            IF (blks(iblk,jblk)%ready) THEN

               IF (allocated(blks(iblk,jblk)%dir  )) deallocate (blks(iblk,jblk)%dir  )
               IF (allocated(blks(iblk,jblk)%upa  )) deallocate (blks(iblk,jblk)%upa  )
               IF (allocated(blks(iblk,jblk)%elv  )) deallocate (blks(iblk,jblk)%elv  )
               IF (allocated(blks(iblk,jblk)%hnd  )) deallocate (blks(iblk,jblk)%hnd  )
               IF (allocated(blks(iblk,jblk)%icat )) deallocate (blks(iblk,jblk)%icat )
               IF (allocated(blks(iblk,jblk)%hunit)) deallocate (blks(iblk,jblk)%hunit)
               IF (allocated(blks(iblk,jblk)%lake )) deallocate (blks(iblk,jblk)%lake )

            ENDIF
         ENDDO
      ENDDO

      thisinfo => allinfo
      DO WHILE (associated(thisinfo))

         IF (allocated (thisinfo%riv_pix)) deallocate (thisinfo%riv_pix)
         IF (allocated (thisinfo%riv_len)) deallocate (thisinfo%riv_len)
         IF (allocated (thisinfo%riv_elv)) deallocate (thisinfo%riv_elv)
         IF (allocated (thisinfo%riv_dep)) deallocate (thisinfo%riv_dep)

         IF (allocated (thisinfo%hru_indx))   deallocate (thisinfo%hru_indx)
         IF (allocated (thisinfo%hru_area))   deallocate (thisinfo%hru_area)
         IF (allocated (thisinfo%hru_next))   deallocate (thisinfo%hru_next)
         IF (allocated (thisinfo%hru_hand))   deallocate (thisinfo%hru_hand)
         IF (allocated (thisinfo%hru_elva))   deallocate (thisinfo%hru_elva)
         IF (allocated (thisinfo%hru_plen))   deallocate (thisinfo%hru_plen)
         IF (allocated (thisinfo%hru_lfac))   deallocate (thisinfo%hru_lfac)

         IF (allocated (thisinfo%bsn_nswe      )) deallocate (thisinfo%bsn_nswe      )
         IF (allocated (thisinfo%bsn_num_hru   )) deallocate (thisinfo%bsn_num_hru   )
         IF (allocated (thisinfo%bsn_downstream)) deallocate (thisinfo%bsn_downstream)
         IF (allocated (thisinfo%bsn_num_nbr   )) deallocate (thisinfo%bsn_num_nbr   )
         IF (allocated (thisinfo%bsn_idx_nbr   )) deallocate (thisinfo%bsn_idx_nbr   )
         IF (allocated (thisinfo%bsn_len_bdr   )) deallocate (thisinfo%bsn_len_bdr   )
         IF (allocated (thisinfo%bsn_elva      )) deallocate (thisinfo%bsn_elva      )

         IF (allocated (thisinfo%lake_id)) deallocate (thisinfo%lake_id)

         thisinfo => thisinfo%next

      ENDDO

   END SUBROUTINE free_memory

   !----------------------------------
   SUBROUTINE append_plist (plist, ilist, i, j, check_exist)
      IMPLICIT NONE

      integer (kind=4), allocatable, intent(inout) :: plist(:,:)
      integer (kind=4), intent(inout) :: ilist
      integer (kind=4), intent(in)    :: i, j
      logical, intent(in) :: check_exist

      integer (kind=4) :: il, nlist
      integer (kind=4), allocatable :: temp(:,:)

      IF (check_exist) THEN
         DO il = 1, ilist
            IF ((plist(1,il) == i) .and. (plist(2,il) == j)) THEN
               RETURN
            ENDIF
         ENDDO
      ENDIF

      ilist = ilist + 1

      IF (.not. allocated(plist)) THEN
         allocate (plist (2,1))
      ENDIF

      nlist = size(plist,2)
      IF (ilist > nlist) THEN
         allocate (temp (2,nlist))
         temp = plist

         deallocate (plist)
         allocate (plist (2,2*nlist))

         plist (:,1:nlist) = temp
      ENDIF

      plist(:,ilist) = (/i, j/)

   END SUBROUTINE append_plist

   !----------------------------------
   SUBROUTINE append_plist3 (plist, ilist, i, j, vlist, val)
      IMPLICIT NONE

      integer (kind=4), allocatable, intent(inout) :: plist(:,:)
      integer (kind=4), intent(inout) :: ilist
      integer (kind=4), intent(in)    :: i, j

      real    (kind=4), allocatable, intent(inout) :: vlist(:)
      real    (kind=4), intent(in) :: val

      integer (kind=4) :: listsize, lz
      integer (kind=4), allocatable :: temp(:,:)
      real    (kind=4), allocatable :: vtmp(:)

      ilist = ilist + 1

      IF (.not. allocated(plist)) THEN
         allocate (plist (2,10))
         allocate (vlist (10))
      ENDIF

      listsize = size(plist,2)
      IF (ilist > listsize) THEN
         lz = max(listsize+10, ceiling(1.2*real(listsize)))
         allocate (temp (2,listsize))
         temp = plist
         deallocate (plist)
         allocate (plist (2,lz))
         plist (:,1:listsize) = temp
         deallocate (temp)
      ENDIF

      listsize = size(vlist)
      IF (ilist > listsize) THEN
         lz = max(listsize+10, ceiling(1.2*real(listsize)))
         allocate (vtmp (listsize))
         vtmp = vlist
         deallocate (vlist)
         allocate (vlist (lz))
         vlist(1:listsize) = vtmp(1:listsize)
         deallocate (vtmp)
      ENDIF

      plist(:,ilist) = (/i, j/)
      vlist(ilist)   = val

   END SUBROUTINE append_plist3

   !---------------------------------------
   integer(kind=4) FUNCTION min_west (i1, i2)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i1, i2

      IF (max(i1,i2) - min(i1,i2) > mglb/2) THEN
         min_west = max(i1,i2)
      ELSE
         min_west = min(i1,i2)
      ENDIF

   END FUNCTION min_west

   !---------------------------------------
   integer(kind=4) FUNCTION max_east (i1, i2)

      IMPLICIT NONE

      integer (kind=4), intent(in) :: i1, i2

      IF (max(i1,i2) - min(i1,i2) > mglb/2) THEN
         max_east = min(i1,i2)
      ELSE
         max_east = max(i1,i2)
      ENDIF

   END FUNCTION max_east

END MODULE hydro_data_mod
