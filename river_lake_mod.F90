MODULE river_lake_mod

CONTAINS

   subroutine get_river_lake (catsize, catsizemin, rivermouthfile)

      use task_mod
      use hydro_data_mod
      USE utils_mod
      implicit none

      real    (kind=4),   intent(in) :: catsize    ! max catchment size
      real    (kind=4),   intent(in) :: catsizemin ! min catchment size
      character(len=256), intent(in), optional :: rivermouthfile

      ! Local Variables
      integer (kind=4) :: nmouth, imouth
      integer (kind=4), allocatable :: rivermouth(:,:)
      real    (kind=4), allocatable :: upa_mouth (:)

      integer (kind=4) :: iproc, scnt, sdsp, idata
      integer (kind=4), allocatable :: rcnt(:), rdsp(:)

      integer (kind=4) :: ncat, icat, icat_dsp, mouthcat(2)
      integer (kind=4) :: nriv, iriv, head, tail
      integer (kind=4), allocatable :: river    (:,:)
      integer (kind=4), allocatable :: mouthlist(:,:)
      real    (kind=4), allocatable :: rivlen(:), rivelv(:)

      real    (kind=4) :: upa_up(8), area, rlen, rlen_max
      integer (kind=4) :: nbranch, ibranch

      integer (kind=4) :: ij_this(2), i, j, i_up, j_up, i_dn, j_dn, idir
      integer (kind=1) :: dir_this
      integer (kind=4) :: iblk, jblk, iblk0, jblk0

      logical :: is_lake
      integer (kind=4) :: lakeid, nplake, iplake
      integer (kind=4), allocatable :: lakelist  (:,:)

      logical :: check_dup_mouth
      integer (kind=4), allocatable :: order(:)
      logical, allocatable :: mouthmask(:)

      integer :: stat

      check_dup_mouth = .not. present(rivermouthfile)

      if (p_is_master) then

         write(*,'(/A/)') 'Step 1 : Finding all rivers and lake outlets ...'

         IF (present(rivermouthfile)) THEN
            
            open (10, status='OLD', file=rivermouthfile, form="FORMATTED")

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
         ELSE
            allocate (rcnt (0:p_ndata))
            allocate (rdsp (0:p_ndata))

            nmouth = 0
            call mpi_gather (nmouth, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, 0, p_comm_data, p_err)

            rdsp(0) = 0
            do iproc = 1, p_ndata
               rdsp(iproc) = rdsp(iproc-1) + rcnt(iproc-1)
            end do

            nmouth = sum(rcnt)
            allocate (rivermouth (2,nmouth))
            allocate (upa_mouth  (nmouth))
            call mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
               rivermouth, rcnt*2, rdsp*2, MPI_INTEGER, 0, p_comm_data, p_err)
            call mpi_gatherv (MPI_IN_PLACE, 0, MPI_REAL4, &
               upa_mouth, rcnt, rdsp, MPI_REAL4, 0, p_comm_data, p_err)

            deallocate (rcnt)
            deallocate (rdsp)
         ENDIF

         call mpi_send (nmouth,     1,      MPI_INTEGER, 1, 1, p_comm_work, p_err)
         call mpi_send (rivermouth, nmouth*2, MPI_INTEGER, 1, 1, p_comm_work, p_err)
         deallocate (rivermouth)

         IF (check_dup_mouth) THEN 
            call mpi_send (upa_mouth, nmouth, MPI_REAL4, 1, 1, p_comm_work, p_err)
            deallocate (upa_mouth)
         ENDIF

         ! river length and elevations.
         call mpi_recv (nrivseg, 1, MPI_INTEGER, 1, 1, p_comm_work, p_stat, p_err)
         allocate (riv_info_len (nrivseg))
         allocate (riv_info_elv (nrivseg))
         call mpi_recv (riv_info_len, nrivseg, MPI_REAL4, 1, 1, p_comm_work, p_stat, p_err)
         call mpi_recv (riv_info_elv, nrivseg, MPI_REAL4, 1, 1, p_comm_work, p_stat, p_err)

         ! gather river pixels.
         call mpi_recv (nrivpix, 1, MPI_INTEGER, 1, 1, p_comm_work, p_stat, p_err)
         allocate (riv_info_pix (3,nrivpix))
         call mpi_recv (riv_info_pix, nrivpix*3, MPI_INTEGER, 1, 1, p_comm_work, p_stat, p_err)

         call excute_data_task (t_exit)

      elseif (p_is_data) then

         IF (.not. present(rivermouthfile)) THEN

            allocate (rivermouth (2,1000000))
            allocate (upa_mouth  (1000000))

            nmouth = 0
            do jblk = 1, mblock
               do iblk = 1, nblock
                  if (bkid(iblk,jblk) == p_iam_glb) then
                     do j = 1, mbox
                        do i = 1, nbox
                           if (blks(iblk,jblk)%icat(i,j) == 0) then
                              ! river mouth and inland depression.
                              IF ((blks(iblk,jblk)%dir(i,j) == 0) .or. (blks(iblk,jblk)%dir(i,j) == -1)) THEN

                                 IF (blks(iblk,jblk)%upa(i,j) >= catsizemin) THEN
                                    call append_plist (rivermouth, nmouth, &
                                       blks(iblk,jblk)%idsp + i, blks(iblk,jblk)%jdsp + j, &
                                       check_exist = .false.)
                                    upa_mouth(nmouth) = blks(iblk,jblk)%upa(i,j)
                                 ENDIF

                              ELSEIF ((blks(iblk,jblk)%idsp+i == inorth)  &
                                 .or. (blks(iblk,jblk)%idsp+i == isouth)  &
                                 .or. (blks(iblk,jblk)%jdsp+j == jwest )  &
                                 .or. (blks(iblk,jblk)%jdsp+j == jeast ))  THEN

                                 IF (blks(iblk,jblk)%upa(i,j) >= catsize) THEN
                                    call append_plist (rivermouth, nmouth, &
                                       blks(iblk,jblk)%idsp + i, blks(iblk,jblk)%jdsp + j, &
                                       check_exist = .false.)
                                    upa_mouth(nmouth) = blks(iblk,jblk)%upa(i,j)
                                 ENDIF

                              ENDIF
                           ENDIF
                        end do
                     end do
                  end if
               end do
            end do

            call mpi_gather  (nmouth, 1, MPI_INTEGER, 0, 0, MPI_INTEGER, 0, p_comm_data, p_err)
            call mpi_gatherv (rivermouth(:,1:nmouth), nmouth*2, MPI_INTEGER, &
               0, 0, 0, MPI_INTEGER, 0, p_comm_data, p_err)
            call mpi_gatherv (upa_mouth(1:nmouth), nmouth, MPI_REAL4, &
               0, 0, 0, MPI_REAL4, 0, p_comm_data, p_err)

            deallocate (rivermouth)
            deallocate (upa_mouth )

         ENDIF

         call data_daemon ()

      end if

      if (p_is_work .and. (p_iam_work == 1)) then

         rlen_max = sqrt(catsize)

         call mpi_recv (nmouth, 1, MPI_INTEGER, 0, 1, p_comm_work, p_stat, p_err)
         
         allocate (rivermouth (2,nmouth))
         call mpi_recv (rivermouth, nmouth*2, MPI_INTEGER, 0, 1, p_comm_work, p_stat, p_err)

         IF (check_dup_mouth) THEN 
            allocate (upa_mouth (nmouth))
            call mpi_recv (upa_mouth, nmouth, MPI_REAL4, 0, 1, p_comm_work, p_stat, p_err)

            allocate (order (nmouth))
            order = (/(imouth, imouth = 1, nmouth)/)

            CALL quicksort (nmouth, upa_mouth, order)

            order = order(nmouth:1:-1)

            rivermouth(1,:) = rivermouth(1,order)
            rivermouth(2,:) = rivermouth(2,order)
            
            deallocate (upa_mouth)
            deallocate (order)
         ENDIF
            
         allocate (river (3,1))
         allocate (mouthlist (2,1000000))
         allocate (lakelist  (2,1000000))

         allocate (mouthmask (nmouth)); mouthmask(:) = .true.

         ncat = 0
         icat = 0
         nriv = 0
         
         do imouth = 1, nmouth

            IF (.not. mouthmask(imouth)) CYCLE

            ! ---------------- next river mouth ---------------------------- 
            mouthcat = rivermouth(:,imouth)
            ! new river segment (case 1): river mouth
            call append_plist (mouthlist, ncat, mouthcat(1), mouthcat(2), check_exist = .false.)

            do while (icat < ncat)

               icat = icat + 1
               mouthcat = mouthlist(:,icat)
               
               is_lake = (get_lake(mouthcat(1),mouthcat(2)) > 0)

               IF (.not. is_lake) THEN

                  write (*,100) mouthcat(1), mouthcat(2), get_upa(mouthcat(1),mouthcat(2))
                  100 format('(S1) River : mouth (',I6,',',I6,'), upa ',F8.0)

                  area = 0
                  rlen = 0

                  ij_this = mouthcat 

                  do while (.true.)

                     IF (get_lake(ij_this(1),ij_this(2)) > 0) THEN
                        ! new river segment (case 2): lake outlet
                           
                        call append_plist (mouthlist, ncat, ij_this(1), ij_this(2), check_exist = .false.)

                        IF (check_dup_mouth) THEN
                           CALL check_duplicate_rivermouth ( &
                              ij_this, rivermouth, nmouth, imouth, mouthmask)
                        ENDIF
                        exit
                     ENDIF

                     upa_up(:) = 0

                     do idir = 1, 8
                        call nextij (ij_this(1), ij_this(2), ishftc(int(-128,1),idir), i_up, j_up)

                        if (((ij_this(1) /= i_up) .or. (ij_this(2) /= j_up)) &
                           .and. within_region(i_up,j_up)) then
                           dir_this = get_dir(i_up,j_up)
                           if (dir_this == ishftc(int(8,1),idir)) then
                              upa_up(idir) = get_upa(i_up,j_up)
                           end if
                        end if
                     end do

                     area = area + sum(upa_up, mask=(upa_up<catsize))
                     nbranch = count(upa_up >= catsize)

                     IF ((area >= catsize) .or. (nbranch /= 1)) THEN
                        ! 3) area reaches threshold; 4) branch point; 5) head water
                     
                        IF ((nbranch == 1) &
                           .and. ((ij_this(1) /= mouthcat(1)) .or. (ij_this(2) /= mouthcat(2)))) then
                           ! new river segment (case 3): reach catchment threshold size.
                              
                           call append_plist (mouthlist, ncat, ij_this(1), ij_this(2), check_exist = .false.)
                              
                           IF (check_dup_mouth) THEN
                              CALL check_duplicate_rivermouth ( &
                                 ij_this, rivermouth, nmouth, imouth, mouthmask)
                           ENDIF
                        ELSE
                           
                           call append_river (river, nriv, ij_this, icat)

                           IF (check_dup_mouth) THEN
                              CALL check_duplicate_rivermouth ( &
                                 ij_this, rivermouth, nmouth, imouth, mouthmask)
                           ENDIF
                          
                           do ibranch = 1, nbranch
                              idir = maxloc(upa_up, dim = 1)
                              call nextij (ij_this(1), ij_this(2), ishftc(int(-128,1),idir), i_up, j_up)

                              ! new river segment (case 4): large branch
                                 
                              call append_plist (mouthlist, ncat, i_up, j_up, check_exist = .false.)
                           
                              IF (check_dup_mouth) THEN
                                 CALL check_duplicate_rivermouth ( &
                                    (/i_up,j_up/), rivermouth, nmouth, imouth, mouthmask)
                              ENDIF

                              upa_up(idir) = 0
                           end do

                           do while ((area >= catsize) .and. any(upa_up > 0))
                              idir = maxloc(upa_up, dim = 1, mask=(upa_up<catsize))
                              call nextij (ij_this(1), ij_this(2), ishftc(int(-128,1),idir), i_up, j_up)

                              ! new river segment (case 5): add small branches 

                              call append_plist (mouthlist, ncat, i_up, j_up, check_exist = .false.)

                              IF (check_dup_mouth) THEN
                                 CALL check_duplicate_rivermouth ( &
                                    (/i_up,j_up/), rivermouth, nmouth, imouth, mouthmask)
                              ENDIF

                              area = area - upa_up(idir)
                              upa_up(idir) = 0
                           end do
                        ENDIF

                        EXIT

                     ELSE ! (area < catsize) .and. (nbranch == 1)

                        ! append to current river
                        call append_river (river, nriv, ij_this, icat)
                                 
                        IF (check_dup_mouth) THEN
                           CALL check_duplicate_rivermouth ( &
                              ij_this, rivermouth, nmouth, imouth, mouthmask)
                        ENDIF

                        idir = maxloc(upa_up, dim = 1)
                        call nextij (ij_this(1), ij_this(2), ishftc(int(-128,1),idir), i_up, j_up)

                        rlen = rlen + dist_between (i_up, j_up, ij_this(1), ij_this(2))

                        if (rlen > rlen_max) then
                           ! new river segment (case 6): too long river
                              
                           call append_plist (mouthlist, ncat, i_up, j_up, check_exist = .true.)

                           IF (check_dup_mouth) THEN
                              CALL check_duplicate_rivermouth ( &
                                 (/i_up,j_up/), rivermouth, nmouth, imouth, mouthmask)
                           ENDIF
                           exit
                        else
                           ij_this = (/i_up,j_up/)
                        end if         

                     ENDIF

                  ENDDO

               ELSE
                  ! lake
                  write (*,101) mouthcat(1), mouthcat(2), get_upa(mouthcat(1),mouthcat(2))
                  101 format('(S1) Lake : outlet (',I6,',',I6,'), upa ',F8.0)

                  lakeid = get_lake(mouthcat(1),mouthcat(2))
                  
                  call append_river (river, nriv, mouthcat, icat)
                                    
                  nplake = 0
                  CALL append_plist (lakelist, nplake, mouthcat(1), mouthcat(2), check_exist = .false.)

                  iplake = 1
                  DO WHILE (iplake <= nplake)
                     do idir = 1, 8
                        call nextij (lakelist(1,iplake), lakelist(2,iplake), &
                           ishftc(int(-128,1),idir), i_up, j_up)

                        if (((lakelist(1,iplake) /= i_up) .or. (lakelist(2,iplake) /= j_up)) &
                           .and. within_region(i_up,j_up)) then
                           if (get_dir(i_up,j_up) == ishftc(int(8,1),idir)) then
                             
                              if (get_lake(i_up,j_up) == lakeid) then
                                 call append_plist (lakelist, nplake, i_up, j_up, check_exist = .false.)
                              ELSEIF (get_upa(i_up,j_up) >= catsize) then
                                 CALL append_plist (mouthlist, ncat, i_up, j_up, check_exist = .false.)
                              ENDIF 
                              
                              IF (get_upa(i_up,j_up) >= catsize) then
                                 IF (check_dup_mouth) THEN
                                    CALL check_duplicate_rivermouth ( &
                                       (/i_up,j_up/), rivermouth, nmouth, imouth, mouthmask)
                                 ENDIF
                              end if

                           end if
                        end if
                     end do

                     iplake = iplake + 1
                  end do                     

               ENDIF
            end do
         end do


         head = 0
         do while (head < nriv)

            tail = head + 1
            iblk = (river(1,tail)-1)/nbox + 1
            jblk = (river(2,tail)-1)/mbox + 1

            head = tail
            do while (head < nriv)
               iblk0 = (river(1,head+1)-1)/nbox + 1
               jblk0 = (river(2,head+1)-1)/mbox + 1
               if ((iblk0 == iblk) .and. (jblk0 == jblk)) then
                  head = head + 1
               else
                  exit
               end if
            end do

            idata = bkid(iblk,jblk)
            call excute_data_task (t_update_river, idata)
            call mpi_send (iblk, 1, MPI_INTEGER4, idata, t_update_river, p_comm_glb, p_err) 
            call mpi_send (jblk, 1, MPI_INTEGER4, idata, t_update_river, p_comm_glb, p_err) 
            call mpi_send (head-tail+1, 1, MPI_INTEGER4, idata, t_update_river, p_comm_glb, p_err) 
            call mpi_send (river(:,tail:head), (head-tail+1)*3, MPI_INTEGER4, &
               idata, t_update_river, p_comm_glb, p_err) 

         end do

         ! river length and elevations.
         allocate (rivlen (ncat))
         allocate (rivelv (ncat))
         icat = 0
         head = 0 
         DO WHILE (icat < ncat)
            tail = head + 1
            head = tail
            DO WHILE (head < nriv)
               if (river(3,head+1) == river(3,head)) then
                  head = head + 1
               else
                  exit
               end if
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

         call mpi_send (ncat,   1,    MPI_INTEGER, 0, 1, p_comm_work, p_err)
         call mpi_send (rivlen, ncat, MPI_REAL4,   0, 1, p_comm_work, p_err)
         call mpi_send (rivelv, ncat, MPI_REAL4,   0, 1, p_comm_work, p_err)


         ! gather river pixels.
         call mpi_send (nriv, 1, MPI_INTEGER, 0, 1, p_comm_work, p_err)
         call mpi_send (river(:,1:nriv), nriv*3, MPI_INTEGER, 0, 1, p_comm_work, p_err)

         IF (allocated(river)     )   deallocate (river)
         IF (allocated(rivermouth))   deallocate (rivermouth)
         IF (allocated(mouthlist ))   deallocate (mouthlist)
         IF (allocated(rivlen)    )   deallocate (rivlen)
         IF (allocated(rivelv)    )   deallocate (rivelv)
         IF (allocated(lakelist  ))   deallocate (lakelist)
         IF (allocated(mouthmask ))   deallocate (mouthmask)

      end if

      call mpi_barrier (p_comm_glb, p_err)

   end subroutine get_river_lake

   ! -- --
   SUBROUTINE check_duplicate_rivermouth (ij, rivmth, nm, im, mmsk)
      
      use hydro_data_mod
      IMPLICIT NONE

      integer, intent(in)    :: ij(2), nm, rivmth(2,nm), im
      logical, intent(inout) :: mmsk(nm)

      ! Local
      integer :: jm

      IF (    (ij(1) == inorth) .or. (ij(1) == isouth)  &
         .or. (ij(2) == jwest ) .or. (ij(2) == jeast ))  THEN
         DO jm= im+1, nm
            IF ((rivmth(1,jm) == ij(1)) .and. (rivmth(2,jm) == ij(2)) ) THEN
               mmsk(jm) = .false.
            ENDIF
         ENDDO
      ENDIF

   END SUBROUTINE check_duplicate_rivermouth

END MODULE river_lake_mod
