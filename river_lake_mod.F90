MODULE river_lake_mod
   
   USE hydro_data_mod

   integer (kind=4) :: nmouth
   integer (kind=4), allocatable :: rivermouth(:,:)
   real    (kind=4), allocatable :: upa_mouth (:)

   integer, SAVE :: ithisblk = -1 
   integer, SAVE :: jthisblk = -1 
   logical :: firstblock = .true.
      
   character(len=256) :: binfo
      
CONTAINS
   
   !----------------------------------------
   logical FUNCTION has_duplicate_mouth (i, j)

      IMPLICIT NONE
      integer, intent(in) :: i, j

      integer :: iglb, jglb, iblk, jblk, iblk_, jblk_, iloc, jloc, inext, jnext
      character(len=256) :: filename
      integer(kind=1), allocatable :: dir  (:,:)  ! flow direction

      iglb = i
      jglb = j
      iblk = (iglb-1)/nbox + 1
      jblk = (jglb-1)/mbox + 1
      iloc = iglb - (iblk-1)*nbox
      jloc = jglb - (jblk-1)*mbox

      CALL get_filename (hydro_dir, iblk, jblk, filename)
      CALL ncio_read_serial (filename, 'dir', dir)

      has_duplicate_mouth = .false.

      DO WHILE ((dir(iloc,jloc) /= 0) .and. (dir(iloc,jloc) /= -1))

         CALL nextij (iglb, jglb, dir(iloc,jloc), inext, jnext)

         IF (on_region_boundary(inext,jnext)) THEN
            has_duplicate_mouth = .true.
            EXIT
         ELSE
            iglb = inext
            jglb = jnext

            iblk_ = (iglb-1)/nbox + 1
            jblk_ = (jglb-1)/mbox + 1
            IF ((iblk_ /= iblk) .or. (jblk_ /= jblk)) THEN
               iblk = iblk_
               jblk = jblk_
               CALL get_filename (hydro_dir, iblk, jblk, filename)
               CALL ncio_read_serial (filename, 'dir', dir)
            ENDIF

            iloc = iglb - (iblk-1)*nbox
            jloc = jglb - (jblk-1)*mbox
         ENDIF
      ENDDO

      deallocate (dir)

   END FUNCTION 


   SUBROUTINE readin_rivermouth (rmfile)

      IMPLICIT NONE

      character(len=*), intent(in) :: rmfile
      ! local
      integer :: stat
      integer (kind=4) :: imouth, mouthcat(2)

      open (10, status='OLD', file=rmfile, form="FORMATTED")

      nmouth = 0
      DO WHILE (.true.)
         read(10,*,iostat=stat) mouthcat
         IF (is_iostat_end(stat)) EXIT
         nmouth = nmouth + 1
      ENDDO

      rewind(10)

      allocate (rivermouth (2,nmouth))
      DO imouth = 1, nmouth
         read(10,*) rivermouth(:,imouth)
      ENDDO
      
      close(10)

   END SUBROUTINE readin_rivermouth

   SUBROUTINE get_next_rivermouth (catsize, catsizemin, end_of_data)

      USE ncio_serial
      IMPLICIT NONE

      real (kind = 4), intent(in) :: catsize, catsizemin
      logical, intent(out) :: end_of_data

      integer(kind=1), allocatable :: dir  (:,:)  ! flow direction
      real   (kind=4), allocatable :: upa  (:,:)  ! upstream drainage area (km^2) 

      type mouth_node
         integer :: i, j
         real    :: upa
         type(mouth_node), pointer :: next
      END type mouth_node

      type(mouth_node), pointer :: allmouth, thismouth

      character(len=256) :: filename, pfmt
      logical :: fexists, mouthfound
      integer :: idsp, jdsp, i, j, m, imouth, jmouth, inext, jnext

      end_of_data = .false.

      DO WHILE (.not. end_of_data)

         IF (firstblock) THEN
            ithisblk = (inorth - 1)/nbox + 1
            jthisblk = (jwest  - 1)/mbox + 1
         ELSE
            idsp = (ithisblk-1)*nbox 
            jdsp = (jthisblk-1)*mbox 
            IF (jdsp > jeast) THEN
               jthisblk = mod(jthisblk,mblock) + 1
            ELSEIF (jdsp+mbox < jeast) THEN
               jthisblk = jthisblk + 1
            ELSEIF (idsp+nbox < isouth) THEN
               ithisblk = ithisblk + 1
               jthisblk = (jwest - 1)/mbox + 1
            ELSE
               end_of_data = .true. 
            ENDIF
         ENDIF

         IF (end_of_data) THEN
            EXIT
         ELSE
            
            idsp = (ithisblk-1)*nbox 
            jdsp = (jthisblk-1)*mbox 

            CALL get_filename (hydro_dir, ithisblk, jthisblk, filename)
            inquire (file=trim(filename), exist=fexists)

            IF (fexists) THEN

               CALL ncio_read_serial (filename, 'dir', dir)
               CALL ncio_read_serial (filename, 'upa', upa)

               nmouth = 0
                  
               write(binfo, '(A,I2,A,I2,A)') '(',ithisblk,',',jthisblk,')'
               write(*,'(/3A/)') 'Step 1 ', trim(binfo), &
                  ': Finding river mouths, inland depressions and outlets on region boundaries ...'

               pfmt = "('(S1) ',A,' mouth found : (',I6,',',I6,'), upa ',F8.0, A)"

               DO j = 1, mbox
                  DO i = 1, nbox
                     imouth = idsp + i
                     jmouth = jdsp + j

                     mouthfound = .false.
                     IF (within_region(imouth, jmouth,.true.)) THEN
                        IF (dir(i,j) == 0) THEN
                           IF (upa(i,j) > catsizemin) THEN
                              mouthfound = .true.
                              write (*,pfmt) trim(binfo), imouth, jmouth, upa(i,j), ' river mouth.'
                           ENDIF
                        ELSEIF (dir(i,j) == -1) THEN
                           mouthfound = .true.
                           write (*,pfmt) trim(binfo), imouth, jmouth, upa(i,j), ' inland depression.'
                        ELSEIF (on_region_boundary(imouth,jmouth) .and. (upa(i,j) >= catsize)) THEN
                           CALL nextij (imouth, jmouth, dir(i,j), inext, jnext)
                           IF (.not. within_region(inext,jnext,.true.)) THEN
                              IF (include_all_upstream) THEN
                                 mouthfound = .not. has_duplicate_mouth (imouth, jmouth)
                              ELSE
                                 mouthfound = .true.
                              ENDIF

                              IF (mouthfound) THEN
                                 write (*,pfmt) trim(binfo), imouth, jmouth, upa(i,j), ' region boundary.'
                              ENDIF
                           ENDIF
                        ENDIF

                        IF (mouthfound) THEN

                           nmouth = nmouth + 1

                           IF (nmouth == 1) THEN
                              allocate(allmouth)
                              thismouth => allmouth
                           ELSE
                              allocate(thismouth%next)
                              thismouth => thismouth%next
                           ENDIF
                           thismouth%next => null()

                           thismouth%i = imouth
                           thismouth%j = jmouth
                           thismouth%upa = upa(i,j)
                        
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO

               IF (nmouth > 0) THEN

                  IF (allocated(rivermouth)) deallocate(rivermouth)
                  IF (allocated(upa_mouth )) deallocate(upa_mouth )

                  allocate (rivermouth (2,nmouth))
                  allocate (upa_mouth    (nmouth))

                  thismouth => allmouth
                  m = 0
                  DO WHILE (associated(thismouth))
                     m = m + 1
                     rivermouth(:,m) = (/thismouth%i, thismouth%j/)
                     upa_mouth   (m) = thismouth%upa

                     allmouth => allmouth%next
                     deallocate (thismouth)
                     thismouth => allmouth
                  ENDDO


                  IF (.not. firstblock) THEN
                     allocate (thisinfo%next)
                     thisinfo%next%icatdsp = thisinfo%icatdsp + thisinfo%ntotalcat
                     thisinfo => thisinfo%next
                     thisinfo%next => null()
                  ENDIF

                  thisinfo%ithisblk = ithisblk
                  thisinfo%jthisblk = jthisblk
         
                  IF (firstblock) firstblock = .false.
               
                  EXIT

               ENDIF
            ENDIF
         ENDIF
      
         IF (firstblock) firstblock = .false.

      ENDDO

      IF (allocated(dir)) deallocate(dir)
      IF (allocated(upa)) deallocate(upa)

   END SUBROUTINE get_next_rivermouth 

   SUBROUTINE get_river_lake (catsize, catsizemin, rivermouthfile, end_of_data)

      USE task_mod
      USE utils_mod
      IMPLICIT NONE

      real    (kind=4),   intent(in) :: catsize    ! max catchment size
      real    (kind=4),   intent(in) :: catsizemin ! min catchment size
      character(len=256), intent(in), optional :: rivermouthfile
      logical, intent(out), optional :: end_of_data

      ! Local Variables
      integer (kind=4) :: imouth, ncat, icat, icat_dsp, mouthcat(2)
      integer (kind=4) :: nriv, iriv, head, tail
      integer (kind=4), allocatable :: river    (:,:)
      integer (kind=4), allocatable :: mouthlist(:,:)
      real    (kind=4), allocatable :: rivlen(:), rivelv(:)

      real    (kind=4) :: upa_up(8), area, rlen, rlen_max
      integer (kind=4) :: nbranch, ibranch

      integer (kind=4) :: ij_this(2), i, j, i_up, j_up, i_dn, j_dn, idir, ip
      integer (kind=1) :: dir_this
      integer (kind=4) :: iblk, jblk, iblk0, jblk0

      logical :: is_lake
      integer (kind=4) :: lakeid, nplake, iplake
      integer (kind=4), allocatable :: lakelist  (:,:)


      IF (p_is_master) THEN
         IF (present(rivermouthfile)) THEN
            CALL readin_rivermouth (rivermouthfile)
            binfo = '(file)'
         ELSE
            CALL get_next_rivermouth (catsize, catsizemin, end_of_data)
            write(binfo, '(A,I2,A,I2,A)') '(',ithisblk,',',jthisblk,')'
         ENDIF
      ENDIF

      IF (present(end_of_data)) THEN
         CALL mpi_bcast (end_of_data, 1, MPI_LOGICAL, p_master_address, p_comm_glb, p_err)
         IF (end_of_data) RETURN
      ENDIF
      
      IF (p_is_master)  write(*,'(/3A/)') 'Step 1 ', trim(binfo), ': Finding all reach and lake outlets ...'


      IF (p_is_master) THEN

         rlen_max = sqrt(catsize)

         allocate (river (3,1))
         allocate (mouthlist (2,1000000))
         allocate (lakelist  (2,1000000))

         ncat = 0
         icat = 0
         nriv = 0
         
         DO imouth = 1, nmouth

            ! ---------------- next river mouth ---------------------------- 
            mouthcat = rivermouth(:,imouth)
            ! new river segment (CASE 1): river mouth
            CALL append_plist (mouthlist, ncat, mouthcat(1), mouthcat(2), check_exist = .false.)

            DO WHILE (icat < ncat)

               icat = icat + 1
               mouthcat = mouthlist(:,icat)
               
               is_lake = (get_lake(mouthcat(1),mouthcat(2)) > 0)

               IF (.not. is_lake) THEN

                  write (*,100) trim(binfo), mouthcat(1), mouthcat(2), get_upa(mouthcat(1),mouthcat(2))
                  100 format('(S1) Reach ',A,': outlet (',I6,',',I6,'), upa ',F8.0)

                  area = 0
                  rlen = 0

                  ij_this = mouthcat 

                  DO WHILE (.true.)

                     IF (get_lake(ij_this(1),ij_this(2)) > 0) THEN
                        ! new river segment (CASE 2): lake outlet
                           
                        CALL append_plist (mouthlist, ncat, ij_this(1), ij_this(2), check_exist = .false.)

                        EXIT
                     ENDIF

                     upa_up(:) = 0

                     DO idir = 1, 8
                        CALL nextij (ij_this(1), ij_this(2), ishftc(int(-128,1),idir), i_up, j_up)

                        IF (((ij_this(1) /= i_up) .or. (ij_this(2) /= j_up)) &
                           .and. within_region(i_up,j_up,.false.)) THEN
                           dir_this = get_dir(i_up,j_up)
                           IF (dir_this == ishftc(int(8,1),idir)) THEN
                              upa_up(idir) = get_upa(i_up,j_up)
                           ENDIF
                        ENDIF
                     ENDDO

                     area = area + sum(upa_up, mask=(upa_up<catsize))
                     nbranch = count(upa_up >= catsize)

                     IF ((area >= catsize) .or. (nbranch /= 1)) THEN
                        ! 3) area reaches threshold; 4) branch point; 5) head water
                     
                        IF ((nbranch == 1) &
                           .and. ((ij_this(1) /= mouthcat(1)) .or. (ij_this(2) /= mouthcat(2)))) THEN
                           ! new river segment (CASE 3): reach catchment threshold size.
                              
                           CALL append_plist (mouthlist, ncat, ij_this(1), ij_this(2), check_exist = .false.)
                              
                        ELSE
                           
                           CALL append_river (river, nriv, ij_this, icat+thisinfo%icatdsp)

                           DO ibranch = 1, nbranch
                              idir = maxloc(upa_up, dim = 1)
                              CALL nextij (ij_this(1), ij_this(2), ishftc(int(-128,1),idir), i_up, j_up)

                              ! new river segment (CASE 4): large branch
                                 
                              CALL append_plist (mouthlist, ncat, i_up, j_up, check_exist = .false.)
                           
                              upa_up(idir) = 0
                           ENDDO

                           DO WHILE ((area >= catsize) .and. any(upa_up > 0))
                              idir = maxloc(upa_up, dim = 1, mask=(upa_up<catsize))
                              CALL nextij (ij_this(1), ij_this(2), ishftc(int(-128,1),idir), i_up, j_up)

                              ! new river segment (CASE 5): add small branches 

                              CALL append_plist (mouthlist, ncat, i_up, j_up, check_exist = .false.)

                              area = area - upa_up(idir)
                              upa_up(idir) = 0
                           ENDDO
                        ENDIF

                        EXIT

                     ELSE ! (area < catsize) .and. (nbranch == 1)

                        ! append to current river
                        CALL append_river (river, nriv, ij_this, icat+thisinfo%icatdsp)
                                 
                        idir = maxloc(upa_up, dim = 1)
                        CALL nextij (ij_this(1), ij_this(2), ishftc(int(-128,1),idir), i_up, j_up)

                        rlen = rlen + dist_between (i_up, j_up, ij_this(1), ij_this(2))

                        IF (rlen > rlen_max) THEN
                           ! new river segment (CASE 6): too long river
                              
                           CALL append_plist (mouthlist, ncat, i_up, j_up, check_exist = .true.)

                           EXIT
                        ELSE
                           ij_this = (/i_up,j_up/)
                        ENDIF         

                     ENDIF

                  ENDDO

               ELSE
                  ! lake
                  write (*,101) trim(binfo), mouthcat(1), mouthcat(2), get_upa(mouthcat(1),mouthcat(2))
                  101 format('(S1) Lake  ',A,': outlet (',I6,',',I6,'), upa ',F8.0)

                  lakeid = get_lake(mouthcat(1),mouthcat(2))
                  
                  CALL append_river (river, nriv, mouthcat, icat+thisinfo%icatdsp)
                                    
                  nplake = 0
                  CALL append_plist (lakelist, nplake, mouthcat(1), mouthcat(2), check_exist = .false.)

                  iplake = 1
                  DO WHILE (iplake <= nplake)
                     DO idir = 1, 8
                        CALL nextij (lakelist(1,iplake), lakelist(2,iplake), &
                           ishftc(int(-128,1),idir), i_up, j_up)

                        IF (((lakelist(1,iplake) /= i_up) .or. (lakelist(2,iplake) /= j_up)) &
                           .and. within_region(i_up,j_up,.false.)) THEN
                           IF (get_dir(i_up,j_up) == ishftc(int(8,1),idir)) THEN
                             
                              IF (get_lake(i_up,j_up) == lakeid) THEN
                                 CALL append_plist (lakelist, nplake, i_up, j_up, check_exist = .false.)
                              ELSEIF (get_upa(i_up,j_up) >= catsize) THEN
                                 CALL append_plist (mouthlist, ncat, i_up, j_up, check_exist = .false.)
                              ENDIF 
                              
                           ENDIF
                        ENDIF
                     ENDDO

                     iplake = iplake + 1
                  ENDDO                     

               ENDIF
            ENDDO
         ENDDO

         DO ip = 1, nriv
            iblk = (river(1,ip)-1)/nbox + 1
            jblk = (river(2,ip)-1)/mbox + 1
         
            IF (.not. blks(iblk,jblk)%ready)  CALL readin_blocks(iblk, jblk)
            
            i = river(1,ip) - blks(iblk,jblk)%idsp
            j = river(2,ip) - blks(iblk,jblk)%jdsp
            IF ((blks(iblk,jblk)%icat (i,j) /= 0) &
               .and. (blks(iblk,jblk)%icat (i,j) /= river(3,ip))) THEN
               write(*,*) 'mismatch while update catnum', &
                  blks(iblk,jblk)%icat(i,j), river(:,ip)
               STOP
            ENDIF
            blks(iblk,jblk)%icat (i,j) = river(3,ip)
         ENDDO

         ! river length and elevations.
         allocate (rivlen (ncat))
         allocate (rivelv (ncat))
         icat = 0
         head = 0 
         DO WHILE (icat < ncat)
            tail = head + 1
            head = tail
            DO WHILE (head < nriv)
               IF (river(3,head+1) == river(3,head)) THEN
                  head = head + 1
               ELSE
                  EXIT
               ENDIF
            ENDDO

            icat = icat + 1
            rivlen(icat) = sqrt(dlon(river(1,tail))*dlat(river(1,tail))/3.1415926)
            DO iriv = tail, head-1
               rivlen(icat) = rivlen(icat) &
                  + dist_between(river(1,iriv), river(2,iriv), river(1,iriv+1), river(2,iriv+1))
            ENDDO
            rivlen(icat) = rivlen(icat) + sqrt(dlon(river(1,head))*dlat(river(1,head))/3.1415926)

            rivelv(icat) = 0
            DO iriv = tail, head
               rivelv(icat) = rivelv(icat) + get_elv(river(1,iriv),river(2,iriv))
            ENDDO
            rivelv(icat) = rivelv(icat) / (head-tail+1)

         ENDDO

         ! river length and elevations.
         thisinfo%nrivseg = ncat
         allocate (thisinfo%riv_len (ncat))
         allocate (thisinfo%riv_elv (ncat))
         thisinfo%riv_len = rivlen
         thisinfo%riv_elv = rivelv

         ! river pixels.
         thisinfo%nrivpix = nriv
         allocate (thisinfo%riv_pix (3,nriv))
         thisinfo%riv_pix = river(:,1:nriv)

         IF (allocated(river)     )   deallocate (river)
         IF (allocated(rivermouth))   deallocate (rivermouth)
         IF (allocated(mouthlist ))   deallocate (mouthlist)
         IF (allocated(rivlen)    )   deallocate (rivlen)
         IF (allocated(rivelv)    )   deallocate (rivelv)
         IF (allocated(lakelist  ))   deallocate (lakelist)

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

   END SUBROUTINE get_river_lake

END MODULE river_lake_mod
