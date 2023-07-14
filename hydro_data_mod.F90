module hydro_data_mod

   implicit none

   integer (kind=4), parameter :: nbox = 6000
   integer (kind=4), parameter :: mbox = 6000
   integer (kind=4), parameter :: nglb = 216000
   integer (kind=4), parameter :: mglb = 432000
   integer (kind=4), parameter :: nblock = 36
   integer (kind=4), parameter :: mblock = 72


   type block_typ

      integer :: idsp, jdsp

      REAL (kind=8), allocatable :: lat(:)
      REAL (kind=8), allocatable :: lon(:)

      ! Flow direction is prepared in 1-byte SIGNED integer (int8) and is defined as:
      !    1: east, 2: southeast, 4: south, 8: southwest, 16: west, 32: northwest, 64: north. 128: northeast
      !    0: river mouth, -1: inland depression, -9: undefined (ocean)
      ! If a flow direction file is opened as UNSIGNED integer, undefined=247 and inland depression=255
      integer(kind=1), allocatable :: dir  (:,:)  ! flow direction
      real   (kind=4), allocatable :: upa  (:,:)  ! upstream drainage area (km^2) 
      real   (kind=4), allocatable :: elv  (:,:)  ! elevation (m)
      real   (kind=4), allocatable :: wth  (:,:)  ! river width (m)

      real   (kind=4), allocatable :: hnd  (:,:)  ! height above nearest drainage (m)

      ! normal catchment with upstream area larger than maximum catchment size (>0); 
      ! flow to ocean(-9); flow to inland depression(-1);
      ! flow to small river mouth (-2); flow outside the region(-3)
      integer(kind=4), allocatable :: icat (:,:)  ! catchment number
      
      integer(kind=4), allocatable :: hunit(:,:)  ! hillslope units

      integer(kind=4), allocatable :: lake (:,:)  ! lake id

   end type block_typ


   type (block_typ) :: blks (nblock,mblock)
   integer (kind=4) :: bkid (nblock,mblock)

   integer (kind=4) :: inorth, isouth, jwest, jeast
   
   real (kind=4) :: dlat(nglb)
   real (kind=4) :: dlon(nglb)
   
   type (block_typ), PRIVATE :: win(2,2)

   ! riv_info_pix(1:2,:) : river pixel coordinate;  
   ! riv_info_pix(3,:)   : catchment id this pixel located in
   integer (kind=4) :: nrivpix, nrivseg
   integer (kind=4), allocatable :: riv_info_pix (:,:)
   real    (kind=4), allocatable :: riv_info_len (:)
   real    (kind=4), allocatable :: riv_info_elv (:)

   ! hydrounit information
   integer (kind=4), allocatable :: hru_info_indx (:,:)    
   real    (kind=4), allocatable :: hru_info_area (:,:)
   real    (kind=4), allocatable :: hru_info_hand (:,:)
   real    (kind=4), allocatable :: hru_info_elva (:,:)
   integer (kind=4), allocatable :: hru_info_next (:,:)    
   real    (kind=4), allocatable :: hru_info_plen (:,:)
   real    (kind=4), allocatable :: hru_info_lfac (:,:)

   ! basin information
   integer (kind=4) :: ntotalcat
   integer (kind=4), allocatable :: bsn_info_bnds     (:,:) 
   integer (kind=4), allocatable :: bsn_info_downstream (:) 
   integer (kind=4), allocatable :: bsn_info_num_nbr    (:)    
   integer (kind=4), allocatable :: bsn_info_idx_nbr  (:,:)    
   real    (kind=4), allocatable :: bsn_info_len_bdr  (:,:)

   ! lake information
   integer (kind=4), allocatable :: lake_info_id (:)    

   ! ------ task id ------
   integer (kind=4), parameter :: t_exit = 0
   integer (kind=4), parameter :: t_inquire_gridinfo = 1
   integer (kind=4), parameter :: t_inquire_dir      = 10
   integer (kind=4), parameter :: t_inquire_upa      = 11
   integer (kind=4), parameter :: t_inquire_elv      = 12
   integer (kind=4), parameter :: t_inquire_icat     = 13
   integer (kind=4), parameter :: t_inquire_wth      = 14
   integer (kind=4), parameter :: t_inquire_hnd      = 15
   integer (kind=4), parameter :: t_inquire_hunit    = 16
   
   integer (kind=4), parameter :: t_inquire_lake     = 17
   
   integer (kind=4), parameter :: t_update_river     = 21
   integer (kind=4), parameter :: t_update_icatch    = 22
   integer (kind=4), parameter :: t_update_hunit     = 23 

contains

   !----------------------------------------
   subroutine get_filename (hydro_dir, iblk, jblk, filename)
      
      implicit none

      character(len=*), intent(in)  :: hydro_dir
      integer (kind=4), intent(in)  :: iblk, jblk
      character(len=*), intent(out) :: filename
      character(len=3)   :: pre1
      character(len=4)   :: pre2

      if (iblk <= 18) then
         write (pre1,'(A1,I2.2)') 'n', (18-iblk)*5     
      else
         write (pre1,'(A1,I2.2)') 's', (iblk-18)*5
      end if

      if (jblk <= 36) then
         write (pre2,'(A1,I3.3)') 'w', (37-jblk)*5
      else
         write (pre2,'(A1,I3.3)') 'e', (jblk-37)*5
      end if

      filename = trim(hydro_dir) // '/' // trim(pre1) // trim(pre2) // '.nc'


   end subroutine get_filename 

   !----------------------------------------
   subroutine get_region (hydro_dir, west, east, north, south)

      use task_mod
      implicit none

      character(len=*), intent(in) :: hydro_dir
      real (kind=4), intent(in) :: west,  east    ! -180 to 180
      real (kind=4), intent(in) :: north, south   ! -90 to 90

      integer (kind=4) :: iblk, jblk

      character(len=256) :: filename
      logical :: fexist

      if (p_is_master) then
         
         if (west == east) then
            jwest = 0
            jeast = 0
         else
            jwest = floor((west + 180.0004167) * 1200) + 1
            jeast = ceiling((east + 180.0004167) * 1200)

            jwest = max(1,min(jwest,mglb))
            jeast = max(1,min(jeast,mglb))
         end if

         inorth = floor((89.9995833 - north) * 1200) + 1
         isouth = ceiling((89.9995833 - south) * 1200)

         inorth = max(1,min(inorth,nglb))
         isouth = max(1,min(isouth,nglb))

         write(*,'(A,4F8.2)') 'Step 0 : Region (n,s,w,e): ', north, south, west, east
         write(*,'(A,4I8)')   'Step 0 : Iregion(n,s,w,e): ', inorth, isouth, jwest, jeast

         do jblk = 1, mblock

            if (jwest > jeast) then
               if ((jwest > jblk*mbox) .and. (jeast < (jblk-1)*mbox+1)) then
                  bkid(:,jblk) = -1
                  cycle
               end if
            elseif (jwest < jeast) then
               if ((jwest > jblk*mbox) .or. (jeast < (jblk-1)*mbox+1)) then
                  bkid(:,jblk) = -1
                  cycle
               end if
            end if

            do iblk = 1, nblock

               if ((inorth <= iblk*nbox) .and. (isouth >= (iblk-1)*nbox+1)) then

                  call get_filename (hydro_dir, iblk, jblk, filename)
                  inquire (file=filename, exist=fexist)

                  if (fexist) then
                     bkid(iblk,jblk) = 0
                  else
                     bkid(iblk,jblk) = -1
                  end if

               else
                  bkid(iblk,jblk) = -1
               end if 
            end do
         end do

      end if

      call mpi_bcast (inorth, 1, MPI_INTEGER, 0, p_comm_glb, p_err) 
      call mpi_bcast (isouth, 1, MPI_INTEGER, 0, p_comm_glb, p_err) 
      call mpi_bcast (jwest,  1, MPI_INTEGER, 0, p_comm_glb, p_err) 
      call mpi_bcast (jeast,  1, MPI_INTEGER, 0, p_comm_glb, p_err) 

      call mpi_bcast (bkid, nblock*mblock, MPI_INTEGER, 0, p_comm_glb, p_err) 


   end subroutine get_region

   !----------------------------------------
   subroutine readin_blocks (hydro_dir, lake_dir)

      use task_mod
      use ncio_serial
      implicit none

      character(len=*), intent(in) :: hydro_dir
      character(len=*), intent(in) :: lake_dir
      integer (kind=4) :: iproc, iblk, jblk, i, j
      integer (kind=4) :: iwin, jwin

      character(len=256) :: filename, lakefile

      iproc = 0
      do jblk = 1, mblock
         do iblk = 1, nblock

            blks(iblk,jblk)%idsp = (iblk-1)*nbox 
            blks(iblk,jblk)%jdsp = (jblk-1)*mbox 

            if (bkid(iblk,jblk) == 0) then

               iproc = mod(iproc, p_ndata) + 1
               bkid(iblk,jblk) = p_data_list(iproc) 

               if (p_is_master) then
                  call get_filename (hydro_dir, iblk, jblk, filename)
                  write(*,*) trim(filename), ' is on data processor ', bkid(iblk,jblk)
               end if

               if (bkid(iblk,jblk) == p_iam_glb) then

                  allocate (blks(iblk,jblk)%dir (nbox,mbox))
                  allocate (blks(iblk,jblk)%upa (nbox,mbox))
                  allocate (blks(iblk,jblk)%elv (nbox,mbox))

                  allocate (blks(iblk,jblk)%hnd  (nbox,mbox))
                  allocate (blks(iblk,jblk)%icat (nbox,mbox))
                  allocate (blks(iblk,jblk)%hunit(nbox,mbox))
                  
                  allocate (blks(iblk,jblk)%lake (nbox,mbox))

                  call get_filename (hydro_dir, iblk, jblk, filename)
                  call ncio_read_serial (filename, 'latitude',  blks(iblk,jblk)%lat)
                  call ncio_read_serial (filename, 'longitude', blks(iblk,jblk)%lon)

                  call ncio_read_serial (filename, 'dir', blks(iblk,jblk)%dir)
                  call ncio_read_serial (filename, 'upa', blks(iblk,jblk)%upa)
                  call ncio_read_serial (filename, 'elv', blks(iblk,jblk)%elv)
                  call ncio_read_serial (filename, 'wth', blks(iblk,jblk)%wth)
                  
                  call get_filename (lake_dir, iblk, jblk, lakefile)
                  call ncio_read_serial (lakefile, 'lake', blks(iblk,jblk)%lake)

                  blks(iblk,jblk)%icat  = -3

                  do i = 1, nbox
                     if ((blks(iblk,jblk)%idsp + i >= inorth) & 
                        .and. (blks(iblk,jblk)%idsp + i <= isouth)) then

                        do j = 1, mbox
                           if ((jwest > jeast) .and. (blks(iblk,jblk)%jdsp + j < jwest) &
                              .and. (blks(iblk,jblk)%jdsp + j > jeast)) then
                              cycle
                           end if

                           if ((jwest < jeast) .and. ((blks(iblk,jblk)%jdsp + j < jwest) &
                              .or. (blks(iblk,jblk)%jdsp + j > jeast))) then
                              cycle
                           end if

                           if (blks(iblk,jblk)%dir(i,j) == -9) then
                              ! ocean
                              blks(iblk,jblk)%icat(i,j) = -9
                           else
                              blks(iblk,jblk)%icat(i,j) = 0
                           end if

                        end do
                     end if
                  end do

                  blks(iblk,jblk)%hnd   = 0
                  blks(iblk,jblk)%hunit = -1
               end if

            end if 
         end do
      end do

      do jwin = 1, 2
         do iwin = 1, 2
            win(iwin,jwin)%idsp = -nbox
            win(iwin,jwin)%jdsp = -mbox
         end do
      end do
      
   end subroutine readin_blocks


   !------------------------------------------------
   subroutine init_area ()

      use MathConstants, only : deg2rad
      implicit none

      real(kind=8), parameter :: re = 6.37122e3 ! kilometer
      real(kind=8) :: lat0, lat1
      integer (kind=4) :: ilat

      do ilat = 1, nglb
         lat1 = 90.0 - 1 / 1200.0 / 2.0 - 180.0 / nglb *  (ilat-1)
         lat0 = lat1 - 1 / 1200.0
         lat1 = deg2rad * min( 90.0, lat1)
         lat0 = deg2rad * max(-90.0, lat0)

         dlat(ilat) = re * (lat1 - lat0)
         dlon(ilat) = re * (sin(lat1)-sin(lat0)) / (lat1-lat0) * deg2rad/1200.0
      end do

   end subroutine init_area

   !-------------------------------------------------
   subroutine nextij (i, j, dir, inext, jnext, is_local)

      implicit none

      integer(kind=4), intent(in) :: i, j
      integer(kind=1), intent(in) :: dir

      integer(kind=4), intent(out)  :: inext, jnext

      LOGICAL, optional, intent(in) :: is_local

      select case (dir)
      case (1)
         inext = i
         jnext = j + 1
      case (2) 
         inext = i + 1
         jnext = j + 1
      case (4)
         inext = i + 1
         jnext = j
      case (8)
         inext = i + 1
         jnext = j - 1
      case (16)
         inext = i
         jnext = j - 1
      case (32)
         inext = i - 1
         jnext = j - 1
      case (64)
         inext = i - 1
         jnext = j
      case (-128)
         inext = i - 1
         jnext = j + 1
      case default
         inext = i
         jnext = j
      end select

      IF (.not. present(is_local)) THEN
         if (jnext == 0)   jnext = mglb
         if (jnext > mglb) jnext = 1
      ELSEIF (.not. is_local) THEN
         if (jnext == 0)   jnext = mglb
         if (jnext > mglb) jnext = 1
      ENDIF

   end subroutine nextij 

   !----------------------------------------------
   logical function is_feasible_step (i, j, i_dn, j_dn) 

      implicit none

      integer (kind=4), intent(in) :: i, j, i_dn, j_dn

      if ((i /= i_dn) .or. (j /= j_dn)) then
         is_feasible_step = .true.
      else
         is_feasible_step = .false.
         return
      end if

      if ((i_dn < inorth) .or. (i_dn > isouth)) then
         is_feasible_step = .false.
         return
      end if

      if (jwest /= 0) then
         if (jwest < jeast) then
            if ((j_dn < jwest) .or. (j_dn > jeast)) then
               is_feasible_step = .false.
               return
            end if
         else
            if ((j_dn < jwest) .and. (j_dn > jeast)) then
               is_feasible_step = .false.
               return
            end if
         end if
      end if

   end function is_feasible_step

   !---------------------------------------------
   subroutine block_iterator ( & 
         imin, imax, jmin, jmax, iblk, jblk, &
         i0, i1, j0, j1, il0, il1, jl0, jl1, &
         iblk1, jblk1, end_of_data)

      implicit none

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

      if (jmax > blks(iblk,jblk)%jdsp) then
         j1 = min(jmax-blks(iblk,jblk)%jdsp, mbox)
      else
         j1 = mbox
      end if

      if (jmin <= blks(iblk,jblk)%jdsp+mbox) then
         j0  = max(jmin-blks(iblk,jblk)%jdsp, 1)
         jl0 = j0 + blks(iblk,jblk)%jdsp - jmin + 1
         jl1 = j1 + blks(iblk,jblk)%jdsp - jmin + 1
      else
         j0  = 1
         jl0 = j0 + blks(iblk,jblk)%jdsp + mglb - jmin + 1 
         jl1 = j1 + blks(iblk,jblk)%jdsp + mglb - jmin + 1 
      end if

      iblk1 = iblk
      end_of_data = .false.
      if (blks(iblk,jblk)%jdsp > jmax) then
         jblk1 = mod(jblk,mblock) + 1
      elseif (blks(iblk,jblk)%jdsp+mbox < jmax) then
         jblk1 = jblk + 1
      elseif (blks(iblk,jblk)%idsp+nbox < imax) then
         iblk1 = iblk1 + 1
         jblk1 = (jmin-1)/mbox + 1
      else
         end_of_data = .true. 
      end if

   end subroutine block_iterator
      

   !---------------------------------------------
   subroutine win_inquire_data (iwin, jwin, iblk, jblk)

      use task_mod
      implicit none

      integer (kind=4), intent(in) :: iwin, jwin
      integer (kind=4), intent(in) :: iblk, jblk

      integer (kind=4) :: idata, bnds(6)

      idata = bkid(iblk,jblk)

      win(iwin,jwin)%idsp = (iblk-1) * nbox
      win(iwin,jwin)%jdsp = (jblk-1) * mbox

      if (idata /= -1) then

         bnds = (/iblk,jblk,1,nbox,1,mbox/)

         call excute_data_task(t_inquire_dir, idata)
         call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_dir, p_comm_glb, p_err) 
         call mpi_recv (win(iwin,jwin)%dir, nbox*mbox, MPI_INTEGER1, idata, t_inquire_dir, &
            p_comm_glb, p_stat, p_err)

         call excute_data_task(t_inquire_upa, idata)
         call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_upa, p_comm_glb, p_err) 
         call mpi_recv (win(iwin,jwin)%upa, nbox*mbox, MPI_REAL4, idata, t_inquire_upa, &
            p_comm_glb, p_stat, p_err)

         call excute_data_task(t_inquire_icat, idata)
         call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_icat, p_comm_glb, p_err) 
         call mpi_recv (win(iwin,jwin)%icat, nbox*mbox, MPI_INTEGER4, idata, t_inquire_icat, &
            p_comm_glb, p_stat, p_err)

         call excute_data_task(t_inquire_elv, idata)
         call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_elv, p_comm_glb, p_err) 
         call mpi_recv (win(iwin,jwin)%elv, nbox*mbox, MPI_REAL4, idata, t_inquire_elv, &
            p_comm_glb, p_stat, p_err)

         call excute_data_task(t_inquire_wth, idata)
         call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_wth, p_comm_glb, p_err) 
         call mpi_recv (win(iwin,jwin)%wth, nbox*mbox, MPI_REAL4, idata, t_inquire_wth, &
            p_comm_glb, p_stat, p_err)

         call excute_data_task(t_inquire_lake, idata)
         call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_lake, p_comm_glb, p_err) 
         call mpi_recv (win(iwin,jwin)%lake, nbox*mbox, MPI_INTEGER4, idata, t_inquire_lake, &
            p_comm_glb, p_stat, p_err)
      else
         win(iwin,jwin)%dir  = -9
         win(iwin,jwin)%upa  = 0
         win(iwin,jwin)%icat = -9
         win(iwin,jwin)%elv  = 0
         win(iwin,jwin)%wth  = 0
         win(iwin,jwin)%lake = 0
      end if

   end subroutine win_inquire_data 

   !---------------------------------------------
   subroutine init_window ()

      USE task_mod
      implicit none

      integer (kind=4) :: i, j

      IF (p_is_work) THEN
         do j = 1, 2
            do i = 1, 2
               win(i,j)%idsp = -nbox
               win(i,j)%jdsp = -mbox

               allocate (win(i,j)%dir  (nbox,mbox))
               allocate (win(i,j)%upa  (nbox,mbox))
               allocate (win(i,j)%icat (nbox,mbox))
               allocate (win(i,j)%elv  (nbox,mbox))
               allocate (win(i,j)%wth  (nbox,mbox))
               allocate (win(i,j)%lake (nbox,mbox))
            end do
         end do
      ENDIF

   end subroutine init_window

   !---------------------------------------------
   SUBROUTINE sync_window ()
      
      IMPLICIT NONE

      integer (kind=4) :: i, j, iblk, jblk

      do j = 1, 2
         do i = 1, 2
            if ((win(i,j)%idsp /= -nbox) .and. (win(i,j)%jdsp /= -mbox)) then
               iblk = win(i,j)%idsp/nbox + 1 
               jblk = win(i,j)%jdsp/mbox + 1
               CALL win_inquire_data(i,j,iblk,jblk)
            ENDIF
         end do
      end do

   END SUBROUTINE sync_window

   !---------------------------------------------
   subroutine free_window ()

      implicit none

      integer (kind=4) :: i, j

      do j = 1, 2
         do i = 1, 2
            IF (allocated(win(i,j)%dir ))  deallocate (win(i,j)%dir )
            IF (allocated(win(i,j)%upa ))  deallocate (win(i,j)%upa )
            IF (allocated(win(i,j)%icat))  deallocate (win(i,j)%icat)
            IF (allocated(win(i,j)%elv ))  deallocate (win(i,j)%elv )
            IF (allocated(win(i,j)%wth ))  deallocate (win(i,j)%wth )
            IF (allocated(win(i,j)%lake))  deallocate (win(i,j)%lake)
         end do
      end do

   end subroutine free_window

   !------------------------------------------------
   subroutine copy_block_in_window (block_from, block_to)

      implicit none
      type (block_typ), intent(in)    :: block_from
      type (block_typ), intent(inout) :: block_to

      block_to%idsp = block_from%idsp
      block_to%jdsp = block_from%jdsp

      if ((block_from%idsp /= -nbox) .and. (block_from%jdsp /= -mbox)) then
         block_to%dir  = block_from%dir
         block_to%upa  = block_from%upa
         block_to%icat = block_from%icat
         block_to%elv  = block_from%elv
         block_to%wth  = block_from%wth
         block_to%lake = block_from%lake
      end if
      
   end subroutine copy_block_in_window 

   !---------------------------------------------
   subroutine get_win_ij (i, j, iwin, jwin)

      implicit none

      integer (kind=4), intent(in)    :: i, j
      integer (kind=4), intent(out)   :: iwin, jwin

      integer (kind=4) :: iblk, jblk
      
      if (i == win(1,1)%idsp) then
         call copy_block_in_window (win(1,1), win(2,1))
         call copy_block_in_window (win(1,2), win(2,2))

         iblk = i/nbox
         jblk = win(1,1)%jdsp/mbox + 1
         call win_inquire_data (1, 1, iblk, jblk)

         win(1,2)%idsp = -nbox
         win(1,2)%jdsp = -mbox

         iwin = 1
      elseif ((i > win(1,1)%idsp) .and. (i <= win(1,1)%idsp+nbox)) then
         iwin = 1
      elseif ((i > win(1,1)%idsp+nbox) .and. (i <= win(1,1)%idsp+2*nbox)) then
         iwin = 2
      elseif (i == win(1,1)%idsp+2*nbox+1) then
         if ((win(2,1)%idsp == -nbox) .and. (win(2,1)%jdsp == -mbox)) then
            iblk = win(1,1)%idsp/nbox + 2
            jblk = win(1,1)%jdsp/mbox + 1
            call win_inquire_data (1, 1, iblk, jblk)
         else
            call copy_block_in_window (win(2,1), win(1,1))
         end if
         call copy_block_in_window (win(2,2), win(1,2))
         
         iwin = 2
         win(2,1)%idsp = -nbox;   win(2,1)%jdsp = -mbox
         win(2,2)%idsp = -nbox;   win(2,2)%jdsp = -mbox
      else
         do jwin = 1, 2
            do iwin = 1, 2
               win(iwin,jwin)%idsp = -nbox
               win(iwin,jwin)%jdsp = -mbox
            ENDDO
         ENDDO

         iblk = (i-1)/nbox + 1
         jblk = (j-1)/mbox + 1
         iwin = 1
         jwin = 1
         call win_inquire_data (iwin, jwin, iblk, jblk)
      end if

      if ((j == win(1,1)%jdsp) .or. ((win(1,1)%jdsp == 0) .and. (j == mglb))) then
         call copy_block_in_window (win(1,1), win(1,2))
         call copy_block_in_window (win(2,1), win(2,2))

         iblk = win(1,1)%idsp/nbox + 1
         jblk = j/mbox
         call win_inquire_data (1, 1, iblk, jblk)

         win(2,1)%idsp = -nbox
         win(2,1)%jdsp = -mbox

         jwin = 1
      elseif ((j > win(1,1)%jdsp) .and. (j <= win(1,1)%jdsp+mbox)) then
         jwin = 1
      elseif (((j > win(1,1)%jdsp+mbox) .and. (j <= win(1,1)%jdsp+2*mbox)) &
            .or. ((win(1,1)%jdsp+mbox == mglb) .and. (j <= mbox))) then
         jwin = 2
      elseif (j == mod(win(1,1)%jdsp+2*mbox+1, mglb)) then
         if ((win(1,2)%idsp == -nbox) .and. (win(1,2)%jdsp == -mbox)) then
            iblk = win(1,1)%idsp/nbox + 1
            jblk = win(1,1)%jdsp/mbox + 2 
            if (jblk > mblock) jblk = jblk - mblock
            call win_inquire_data (1, 1, iblk, jblk)
         else
            call copy_block_in_window (win(1,2), win(1,1))
         end if
         call copy_block_in_window (win(2,2), win(2,1))
         
         jwin = 2
         win(1,2)%idsp = -nbox;   win(1,2)%jdsp = -mbox
         win(2,2)%idsp = -nbox;   win(2,2)%jdsp = -mbox
      else
         do jwin = 1, 2
            do iwin = 1, 2
               win(iwin,jwin)%idsp = -nbox
               win(iwin,jwin)%jdsp = -mbox
            ENDDO
         ENDDO

         iblk = (i-1)/nbox + 1
         jblk = (j-1)/mbox + 1
         iwin = 1
         jwin = 1
         call win_inquire_data (iwin, jwin, iblk, jblk)
      end if

      if ((win(iwin,jwin)%idsp == -nbox) .and. (win(iwin,jwin)%jdsp == -mbox)) then
         iblk = win(1,1)%idsp/nbox + iwin
         jblk = win(1,1)%jdsp/mbox + jwin
         if (jblk > mblock) then
            jblk = jblk - mblock
         end if

         call win_inquire_data (iwin, jwin, iblk, jblk)
      end if

   end subroutine get_win_ij


   !----------------------------------------------
   integer(kind=1) function get_dir (i, j)

      implicit none

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iwin, jwin

      call get_win_ij (i, j, iwin, jwin)

      get_dir = win(iwin,jwin)%dir (i-win(iwin,jwin)%idsp, j-win(iwin,jwin)%jdsp)

   end function get_dir 
      

   !----------------------------------------------
   real(kind=4) function get_upa (i, j)

      implicit none

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iwin, jwin

      call get_win_ij (i, j, iwin, jwin)

      get_upa = win(iwin,jwin)%upa (i-win(iwin,jwin)%idsp, j-win(iwin,jwin)%jdsp)

   end function get_upa

   !----------------------------------------------
   real(kind=4) function get_elv (i, j)

      implicit none

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iwin, jwin

      call get_win_ij (i, j, iwin, jwin)

      get_elv = win(iwin,jwin)%elv (i-win(iwin,jwin)%idsp, j-win(iwin,jwin)%jdsp)

   end function get_elv

   !----------------------------------------------
   integer (kind=4) function get_icat (i, j)

      implicit none

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iwin, jwin

      call get_win_ij (i, j, iwin, jwin)
      
      get_icat = win(iwin,jwin)%icat (i-win(iwin,jwin)%idsp, j-win(iwin,jwin)%jdsp)

   end function get_icat

   !----------------------------------------------
   real(kind=4) function get_wth (i, j)

      implicit none

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iwin, jwin

      call get_win_ij (i, j, iwin, jwin)

      get_wth = win(iwin,jwin)%wth (i-win(iwin,jwin)%idsp, j-win(iwin,jwin)%jdsp)

   end function get_wth

   !----------------------------------------------
   integer(kind=4) function get_lake (i, j)

      implicit none

      integer (kind=4), intent(in) :: i, j

      integer (kind=4) :: iwin, jwin

      call get_win_ij (i, j, iwin, jwin)

      get_lake = win(iwin,jwin)%lake (i-win(iwin,jwin)%idsp, j-win(iwin,jwin)%jdsp)

   end function get_lake

   !----------------------------------------------
   real(kind=4) function get_area (i, j)

      implicit none

      integer (kind=4), intent(in) :: i, j

      get_area = dlon(i) * dlat(i)

   end function get_area
   
   !----------------------------------------------
   real(kind=4) function dist_between (i0, j0, i1, j1)

      implicit none

      integer (kind=4), intent(in) :: i0, j0, i1, j1

      if (i0 == i1) then
         dist_between = dlon(i0)
      elseif (j0 == j1) then
         dist_between = 0.5*(dlat(i0) + dlat(i1))
      else
         dist_between = sqrt((0.5*(dlat(i0)+dlat(i1)))**2 + (0.5*(dlon(i0)+dlon(i1)))**2)
      end if

   end function dist_between

   !------------------------------------
   subroutine aggregate_data (imin, imax, jmin, jmax, &
         np, mp, longitude, latitude, icat, dir, hnd, elv, hunit)
      
      use task_mod
      implicit none
   
      integer (kind=4), intent(in) :: imin, imax, jmin, jmax
         
      integer (kind=4), intent(in) :: np, mp
      REAL    (kind=8), intent(inout) :: longitude(mp)
      REAL    (kind=8), intent(inout) :: latitude (np)

      integer (kind=4), intent(inout), optional :: icat  (np,mp)    
      integer (kind=1), intent(inout), optional :: dir   (np,mp)    
      real    (kind=4), intent(inout), optional :: hnd   (np,mp)
      real    (kind=4), intent(inout), optional :: elv   (np,mp)
      integer (kind=4), intent(inout), optional :: hunit (np,mp)    
   
      integer (kind=4) :: idata, dsize
      integer (kind=4) :: iblk, jblk, iblk1, jblk1
      integer (kind=4) :: i0, i1, j0, j1, il0, il1, jl0, jl1
      integer (kind=4) :: bnds(6)
      logical :: end_of_data
      
      iblk = (imin-1)/nbox + 1
      jblk = (jmin-1)/mbox + 1
      do while (.true.)

         call block_iterator (imin, imax, jmin, jmax, iblk, jblk, &
            i0, i1, j0, j1, il0, il1, jl0, jl1, &
            iblk1, jblk1, end_of_data)

         idata = bkid(iblk,jblk)
         if (idata == -1) then
            ! ocean
            if (present(icat)) then
               icat(il0:il1,jl0:jl1) = -9
            end if
         else

            bnds  = (/iblk,jblk,i0,i1,j0,j1/)
            dsize = (j1-j0+1)*(i1-i0+1)
         
            call excute_data_task (t_inquire_gridinfo, idata)
            call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_gridinfo, p_comm_glb, p_err) 
            call mpi_recv (longitude(jl0:jl1), jl1-jl0+1, MPI_DOUBLE, idata, &
               t_inquire_gridinfo, p_comm_glb, p_stat, p_err)
            call mpi_recv (latitude (il0:il1), il1-il0+1, MPI_DOUBLE, idata, &
               t_inquire_gridinfo, p_comm_glb, p_stat, p_err)

            if (present(icat)) then
               call excute_data_task (t_inquire_icat, idata)
               call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_icat, p_comm_glb, p_err) 
               call mpi_recv (icat(il0:il1,jl0:jl1), dsize, MPI_INTEGER, idata, &
                  t_inquire_icat, p_comm_glb, p_stat, p_err)
            end if

            if (present(dir)) then
               call excute_data_task (t_inquire_dir, idata)
               call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_dir, p_comm_glb, p_err) 
               call mpi_recv (dir(il0:il1,jl0:jl1), dsize, MPI_INTEGER1, idata, &
                  t_inquire_dir, p_comm_glb, p_stat, p_err)
            end if

            if (present(hnd)) then
               call excute_data_task (t_inquire_hnd, idata)
               call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_hnd, p_comm_glb, p_err) 
               call mpi_recv (hnd(il0:il1,jl0:jl1), dsize, MPI_REAL4, idata, &
                  t_inquire_hnd, p_comm_glb, p_stat, p_err)
            end if

            if (present(elv)) then
               call excute_data_task (t_inquire_elv, idata)
               call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_elv, p_comm_glb, p_err) 
               call mpi_recv (elv(il0:il1,jl0:jl1), dsize, MPI_REAL4, idata, &
                  t_inquire_elv, p_comm_glb, p_stat, p_err)
            end if

            if (present(hunit)) then
               call excute_data_task (t_inquire_hunit, idata)
               call mpi_send (bnds, 6, MPI_INTEGER, idata, t_inquire_hunit, p_comm_glb, p_err) 
               call mpi_recv (hunit(il0:il1,jl0:jl1), dsize, MPI_INTEGER, idata, &
                  t_inquire_hunit, p_comm_glb, p_stat, p_err)
            end if

         end if


         if (end_of_data) then
            exit
         else
            iblk = iblk1; jblk = jblk1
         end if
      end do

   end subroutine aggregate_data

   !------------------------------------
   subroutine excute_data_task (t_id, idata)

      use task_mod
      implicit none

      integer (kind=4), intent(in) :: t_id
      integer (kind=4), intent(in), optional :: idata
      integer (kind=4) :: sbuf(2) 
      integer (kind=4) :: iproc

      sbuf = (/t_id, p_iam_glb/)

      if (present(idata)) then
         call mpi_send (sbuf, 2, MPI_INTEGER, idata, t_id, p_comm_glb, p_err) 
      else
         do iproc = 1, p_ndata
            call mpi_send (sbuf, 2, MPI_INTEGER, p_data_list(iproc), t_id, p_comm_glb, p_err) 
         end do
      end if

   end subroutine excute_data_task

   !------------------------------------
   subroutine data_daemon ()

      use task_mod
      implicit none

      integer (kind=4) :: rbuf(2), isrc, t_id
   
      integer (kind=4) :: iblk, jblk
      integer (kind=4) :: i, j, ii, jj 
      
      integer (kind=4) :: nnode, inode
      integer (kind=4), allocatable :: route(:,:)
      integer (kind=4) :: icatch
      real    (kind=4), allocatable :: elvdata(:)
      
      integer (kind=4) :: nriv, iriv
      integer (kind=4), allocatable :: river(:,:)

      integer (kind=4) :: i0, i1, j0, j1
      integer (kind=4) :: bnds(6)
      integer (kind=4) :: dsize

      integer (kind=4), allocatable :: hunit_recv(:,:)


      do while (.true.)

         call mpi_recv (rbuf, 2, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, p_comm_glb, p_stat, p_err)
         t_id = rbuf(1)
         isrc = rbuf(2)

         select case (t_id)

         case (t_inquire_gridinfo)

            call mpi_recv (bnds, 6, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)

            iblk = bnds(1);  jblk = bnds(2)
            i0 = bnds(3); i1 = bnds(4); j0 = bnds(5); j1 = bnds(6)
            call mpi_send (blks(iblk,jblk)%lon(j0:j1), j1-j0+1, MPI_DOUBLE, isrc, &
               t_id, p_comm_glb, p_err) 
            call mpi_send (blks(iblk,jblk)%lat(i0:i1), i1-i0+1, MPI_DOUBLE, isrc, &
               t_id, p_comm_glb, p_err) 

         case (t_inquire_dir)

            call mpi_recv (bnds, 6, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)

            iblk = bnds(1);  jblk = bnds(2)
            i0 = bnds(3); i1 = bnds(4); j0 = bnds(5); j1 = bnds(6)
            dsize = (j1-j0+1)*(i1-i0+1)
            call mpi_send (blks(iblk,jblk)%dir(i0:i1,j0:j1), dsize, MPI_INTEGER1, isrc, &
               t_id, p_comm_glb, p_err) 

         case (t_inquire_upa)

            call mpi_recv (bnds, 6, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)

            iblk = bnds(1);  jblk = bnds(2)
            i0 = bnds(3); i1 = bnds(4); j0 = bnds(5); j1 = bnds(6)
            dsize = (j1-j0+1)*(i1-i0+1)
            call mpi_send (blks(iblk,jblk)%upa(i0:i1,j0:j1), dsize, MPI_REAL4, isrc, &
               t_id, p_comm_glb, p_err) 

         case (t_inquire_elv)

            call mpi_recv (bnds, 6, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)

            iblk = bnds(1);  jblk = bnds(2)
            i0 = bnds(3); i1 = bnds(4); j0 = bnds(5); j1 = bnds(6)
            dsize = (j1-j0+1)*(i1-i0+1)
            call mpi_send (blks(iblk,jblk)%elv(i0:i1,j0:j1), dsize, MPI_REAL4, isrc, &
               t_id, p_comm_glb, p_err) 

         case (t_inquire_icat)

            call mpi_recv (bnds, 6, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)

            iblk = bnds(1);  jblk = bnds(2)
            i0 = bnds(3); i1 = bnds(4); j0 = bnds(5); j1 = bnds(6)
            dsize = (j1-j0+1)*(i1-i0+1)
            call mpi_send (blks(iblk,jblk)%icat(i0:i1,j0:j1), dsize, MPI_INTEGER, isrc, &
               t_id, p_comm_glb, p_err) 

         case (t_inquire_wth)

            call mpi_recv (bnds, 6, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)

            iblk = bnds(1);  jblk = bnds(2)
            i0 = bnds(3); i1 = bnds(4); j0 = bnds(5); j1 = bnds(6)
            dsize = (j1-j0+1)*(i1-i0+1)
            call mpi_send (blks(iblk,jblk)%wth(i0:i1,j0:j1), dsize, MPI_REAL4, isrc, &
               t_id, p_comm_glb, p_err) 

         case (t_inquire_hnd)

            call mpi_recv (bnds, 6, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)

            iblk = bnds(1);  jblk = bnds(2)
            i0 = bnds(3); i1 = bnds(4); j0 = bnds(5); j1 = bnds(6)
            dsize = (j1-j0+1)*(i1-i0+1)
            call mpi_send (blks(iblk,jblk)%hnd(i0:i1,j0:j1), dsize, MPI_REAL4, isrc, &
               t_id, p_comm_glb, p_err) 

         case (t_inquire_hunit)

            call mpi_recv (bnds, 6, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)

            iblk = bnds(1);  jblk = bnds(2)
            i0 = bnds(3); i1 = bnds(4); j0 = bnds(5); j1 = bnds(6)
            dsize = (j1-j0+1)*(i1-i0+1)
            call mpi_send (blks(iblk,jblk)%hunit(i0:i1,j0:j1), dsize, MPI_INTEGER, isrc, &
               t_id, p_comm_glb, p_err) 

         case (t_inquire_lake)

            call mpi_recv (bnds, 6, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)

            iblk = bnds(1);  jblk = bnds(2)
            i0 = bnds(3); i1 = bnds(4); j0 = bnds(5); j1 = bnds(6)
            dsize = (j1-j0+1)*(i1-i0+1)
            call mpi_send (blks(iblk,jblk)%lake(i0:i1,j0:j1), dsize, MPI_INTEGER, isrc, &
               t_id, p_comm_glb, p_err) 

         case (t_update_river)

            call mpi_recv (iblk, 1, MPI_INTEGER4, isrc, t_id, p_comm_glb, p_stat, p_err)
            call mpi_recv (jblk, 1, MPI_INTEGER4, isrc, t_id, p_comm_glb, p_stat, p_err)
            call mpi_recv (nriv, 1, MPI_INTEGER4, isrc, t_id, p_comm_glb, p_stat, p_err)
            allocate (river (3,nriv))
            call mpi_recv (river, nriv*3, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)

            do iriv = 1, nriv
               i = river(1,iriv) - blks(iblk,jblk)%idsp
               j = river(2,iriv) - blks(iblk,jblk)%jdsp
               blks(iblk,jblk)%icat(i,j) = river(3,iriv)
            end do

            deallocate (river)

         case (t_update_icatch)

            call mpi_recv (nnode,  1, MPI_INTEGER4, isrc, t_id, p_comm_glb, p_stat, p_err)
            call mpi_recv (iblk,   1, MPI_INTEGER4, isrc, t_id, p_comm_glb, p_stat, p_err)
            call mpi_recv (jblk,   1, MPI_INTEGER4, isrc, t_id, p_comm_glb, p_stat, p_err)
            call mpi_recv (icatch, 1, MPI_INTEGER4, isrc, t_id, p_comm_glb, p_stat, p_err)

            allocate (route (2,nnode))
            call mpi_recv (route, nnode*2, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)
            
            allocate (elvdata (nnode))
            call mpi_recv (elvdata, nnode, MPI_REAL4, isrc, t_id, p_comm_glb, p_stat, p_err)

            do inode = 1, nnode
               i = route(1,inode) - blks(iblk,jblk)%idsp
               j = route(2,inode) - blks(iblk,jblk)%jdsp
               blks(iblk,jblk)%icat (i,j) = icatch
               blks(iblk,jblk)%hnd  (i,j) = elvdata(inode)
            end do

            deallocate (route)
            deallocate (elvdata)

         case (t_update_hunit)

            call mpi_recv (icatch, 1, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)
            call mpi_recv (bnds,   6, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err)

            iblk = bnds(1);  jblk = bnds(2)
            i0 = bnds(3); i1 = bnds(4); j0 = bnds(5); j1 = bnds(6)
            dsize = (j1-j0+1)*(i1-i0+1)
            
            allocate (hunit_recv (i1-i0+1,j1-j0+1))
            call mpi_recv (hunit_recv, dsize, MPI_INTEGER, isrc, t_id, p_comm_glb, p_stat, p_err) 

            do i = i0, i1
               do j = j0, j1
                  if (blks(iblk,jblk)%icat(i,j) == icatch) then
                     blks(iblk,jblk)%hunit(i,j) = hunit_recv(i-i0+1,j-j0+1)
                  end if
               end do
            end do

            deallocate (hunit_recv)

         case (t_exit) 

            exit

         end select 

      end do

   end subroutine data_daemon

   !-------------------------------------------------------
   subroutine append_river (river, nriv, ij0, icatch)

      implicit none
      
      integer (kind=4), intent(inout), allocatable :: river(:,:)
      integer (kind=4), intent(inout) :: nriv
      integer (kind=4), intent(in) :: ij0(2), icatch

      integer (kind=4) :: nsize
      integer (kind=4), allocatable :: temp(:,:)
         
      nriv = nriv + 1
      
      nsize = size(river,2)
      if (nriv > nsize) then
         allocate (temp(3,nsize))
         temp = river
         
         deallocate (river)
         allocate   (river(3,nsize+100000))
         river(:,1:nsize) = temp

         deallocate (temp)
      end if

      river (1:2,nriv) = ij0
      river   (3,nriv) = icatch

   end subroutine append_river
   
   !--------------------------------------------------------------
   subroutine free_memory ()

      use task_mod
      implicit none

      integer (kind=4) :: iblk, jblk
                  

      if (p_is_data) then
         do iblk = 1, nblock
            do jblk = 1, mblock
               if (bkid(iblk,jblk) == p_iam_glb) then

                  deallocate (blks(iblk,jblk)%dir)
                  deallocate (blks(iblk,jblk)%upa)
                  deallocate (blks(iblk,jblk)%elv)
                  deallocate (blks(iblk,jblk)%hnd  )
                  deallocate (blks(iblk,jblk)%icat )
                  deallocate (blks(iblk,jblk)%hunit)
                  deallocate (blks(iblk,jblk)%lake )

               end if
            end do
         end do
      end if

      if (allocated (riv_info_pix)) deallocate (riv_info_pix)
      IF (allocated (riv_info_len)) deallocate (riv_info_len)
      IF (allocated (riv_info_elv)) deallocate (riv_info_elv)
      
      IF (allocated (hru_info_indx))   deallocate (hru_info_indx)
      IF (allocated (hru_info_area))   deallocate (hru_info_area)
      IF (allocated (hru_info_next))   deallocate (hru_info_next)
      IF (allocated (hru_info_hand))   deallocate (hru_info_hand)
      IF (allocated (hru_info_elva))   deallocate (hru_info_elva)
      IF (allocated (hru_info_plen))   deallocate (hru_info_plen)
      IF (allocated (hru_info_lfac))   deallocate (hru_info_lfac)

      IF (allocated (bsn_info_bnds      )) deallocate (bsn_info_bnds      )
      IF (allocated (bsn_info_downstream)) deallocate (bsn_info_downstream)
      IF (allocated (bsn_info_num_nbr   )) deallocate (bsn_info_num_nbr   )
      IF (allocated (bsn_info_idx_nbr   )) deallocate (bsn_info_idx_nbr   )
      IF (allocated (bsn_info_len_bdr   )) deallocate (bsn_info_len_bdr   )
      
      IF (allocated (lake_info_id)) deallocate (lake_info_id)

      call free_window ()

   end subroutine free_memory
   
   !----------------------------------
   subroutine append_plist (plist, ilist, i, j)
      implicit none

      integer (kind=4), allocatable, intent(inout) :: plist(:,:)
      integer (kind=4), intent(inout) :: ilist
      integer (kind=4), intent(in)    :: i, j

      integer (kind=4) :: il, nlist
      integer (kind=4), allocatable :: temp(:,:)

      ! DO il = 1, ilist
      !    IF ((plist(1,il) == i) .and. (plist(2,il) == j)) THEN
      !       RETURN
      !    ENDIF
      ! ENDDO

      ilist = ilist + 1

      IF (.not. allocated(plist)) THEN
         allocate (plist (2,1))
      ENDIF

      nlist = size(plist,2)
      if (ilist > nlist) then
         allocate (temp (2,nlist))
         temp = plist

         deallocate (plist)
         allocate (plist (2,2*nlist))

         plist (:,1:nlist) = temp
      end if

      plist(:,ilist) = (/i, j/)

   end subroutine append_plist

   !----------------------------------
   subroutine append_plist3 (plist, ilist, i, j, vlist, val)
      implicit none

      integer (kind=4), allocatable, intent(inout) :: plist(:,:)
      integer (kind=4), intent(inout) :: ilist
      integer (kind=4), intent(in)    :: i, j
      
      real    (kind=4), allocatable, intent(inout) :: vlist(:)
      real    (kind=4), intent(in) :: val

      integer (kind=4) :: listsize
      integer (kind=4), allocatable :: temp(:,:)
      real    (kind=4), allocatable :: vtmp(:)

      ilist = ilist + 1

      IF (.not. allocated(plist)) THEN
         allocate (plist (2,10))
         allocate (vlist (10))
      ENDIF

      listsize = size(plist,2)
      if (ilist > listsize) then
         allocate (temp (2,listsize))
         temp = plist
         deallocate (plist)
         allocate (plist (2,ceiling(1.2*listsize)))
         plist (:,1:listsize) = temp
         
         allocate (vtmp (listsize))
         vtmp = vlist
         deallocate (vlist)
         allocate (vlist (ceiling(1.2*listsize)))
         vlist(1:listsize) = vtmp(1:listsize)

         deallocate (temp)
         deallocate (vtmp)
      end if

      plist(:,ilist) = (/i, j/)
      vlist(ilist)   = val

   end subroutine append_plist3

   !---------------------------------------
   integer(kind=4) function min_west (i1, i2)

      implicit none

      integer (kind=4), intent(in) :: i1, i2

      if (max(i1,i2) - min(i1,i2) > mglb/2) then
         min_west = max(i1,i2)
      else
         min_west = min(i1,i2)
      end if

   end function min_west

   !---------------------------------------
   integer(kind=4) function max_east (i1, i2)

      implicit none

      integer (kind=4), intent(in) :: i1, i2

      if (max(i1,i2) - min(i1,i2) > mglb/2) then
         max_east = min(i1,i2)
      else
         max_east = max(i1,i2)
      end if

   end function max_east

end module hydro_data_mod
