MODULE neighbour_mod

CONTAINS

   SUBROUTINE get_basin_neighbour (ntotalall, nnbmax)

      USE task_mod
      USE hydro_data_mod

      IMPLICIT NONE

      integer (kind=4), intent(in)  :: ntotalall
      integer (kind=4), intent(out) :: nnbmax

      ! Local Variables
      integer (kind=4), allocatable :: catch (:,:)

      integer (kind=4) :: numblocks, ntotalcat, iblk, jblk
      integer (kind=4) :: catnum, imin, imax, jmin, jmax, icat, jcat, np, mp
      integer (kind=4) :: ithis, jthis, inext, jnext, igthis, jgthis, ignext, jgnext
      integer (kind=4) :: iwork, zero, ndone, mesg(5)

      integer (kind = 4) :: nnb, inb
      integer (kind = 4), allocatable :: nbindex  (:)
      real    (kind = 4), allocatable :: lenborder(:)
      real    (kind = 4) :: blen

      type basin_neighbour_type
         integer :: nnb
         integer, allocatable :: nbr_index(:)
         real*4 , allocatable :: lenborder(:)
      END type basin_neighbour_type

      type(basin_neighbour_type), allocatable :: bsn_nbr(:)


      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN

         write(*,'(/A/)') 'Step 4 : Get neighbour information ...'

         numblocks = 0
         thisinfo => allinfo
         DO WHILE (associated(thisinfo))
            numblocks = numblocks + 1
            thisinfo => thisinfo%next
         ENDDO

         CALL mpi_bcast (numblocks, 1, MPI_INTEGER, p_master_address, p_comm_glb, p_err)

         nnbmax = 0

         thisinfo => allinfo
         DO WHILE (.true.)

            ntotalcat = thisinfo%ntotalcat

            allocate (bsn_nbr (ntotalcat))

            icat = 1
            ndone = 0
            DO WHILE (.true.)

               CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, MPI_ANY_SOURCE, 0, p_comm_glb, p_stat, p_err)

               iwork  = mesg(1)
               catnum = mesg(2)
               IF (catnum > 0) THEN

                  jcat = catnum - thisinfo%icatdsp

                  CALL mpi_recv (nnb, 1, MPI_INTEGER, iwork, 2, p_comm_glb, p_stat, p_err)
                  bsn_nbr(jcat)%nnb = nnb
                  IF (nnb > 0) THEN
                     allocate (bsn_nbr(jcat)%nbr_index(nnb))
                     allocate (bsn_nbr(jcat)%lenborder(nnb))
                     CALL mpi_recv (bsn_nbr(jcat)%nbr_index, nnb, MPI_INTEGER, iwork, 2, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv (bsn_nbr(jcat)%lenborder, nnb, MPI_REAL4  , iwork, 2, p_comm_glb, p_stat, p_err)
                  ENDIF

                  nnbmax = max(nnb, nnbmax)

               ENDIF

               IF (icat <= ntotalcat) THEN

                  catnum = thisinfo%icatdsp+icat

                  CALL mpi_send (catnum, 1, MPI_INTEGER, iwork, 1, p_comm_glb, p_err)

                  imin = thisinfo%bsn_nswe(1,icat)
                  imax = thisinfo%bsn_nswe(2,icat)
                  jmin = thisinfo%bsn_nswe(3,icat)
                  jmax = thisinfo%bsn_nswe(4,icat)

                  CALL enlarge_nswe (imin, imax, jmin, jmax)

                  np = imax - imin + 1
                  mp = jmax - jmin + 1
                  IF (mp < 0) mp = mp + mglb

                  allocate (catch (np,mp))

                  CALL aggregate_catch (imin, imax, jmin, jmax, np, mp, catch)

                  mesg(1:4) = (/imin, jmin, np, mp/)
                  CALL mpi_send (mesg(1:4), 4, MPI_INTEGER,  iwork, 1, p_comm_glb, p_err)
                  CALL mpi_send (catch, np*mp, MPI_INTEGER,  iwork, 1, p_comm_glb, p_err)

                  deallocate (catch)

                  icat = icat + 1
               ELSE
                  zero = 0
                  CALL mpi_send (zero, 1, MPI_INTEGER, iwork, 1, p_comm_glb, p_err)
                  ndone = ndone + 1
               ENDIF

               IF (ndone == p_nwork) EXIT
            ENDDO

            allocate (thisinfo%bsn_num_nbr (ntotalcat))
            allocate (thisinfo%bsn_idx_nbr (nnbmax,ntotalcat));   thisinfo%bsn_idx_nbr(:,:) = -1
            allocate (thisinfo%bsn_len_bdr (nnbmax,ntotalcat));   thisinfo%bsn_len_bdr(:,:) = 0

            DO icat = 1, ntotalcat
               thisinfo%bsn_num_nbr(icat) = bsn_nbr(icat)%nnb
               IF (bsn_nbr(icat)%nnb > 0) THEN
                  thisinfo%bsn_idx_nbr(1:bsn_nbr(icat)%nnb,icat) = bsn_nbr(icat)%nbr_index
                  thisinfo%bsn_len_bdr(1:bsn_nbr(icat)%nnb,icat) = bsn_nbr(icat)%lenborder
               ENDIF
            ENDDO

            DO icat = 1, ntotalcat
               IF (bsn_nbr(icat)%nnb > 0) THEN
                  deallocate (bsn_nbr(icat)%nbr_index)
                  deallocate (bsn_nbr(icat)%lenborder)
               ENDIF
            ENDDO
            deallocate (bsn_nbr)

            IF (trim(def_storage_type) == 'block') THEN
               DO jblk = 1, mblock
                  DO iblk = 1, nblock
                     IF (allocated(blks(iblk,jblk)%icat)) THEN
                        deallocate (blks(iblk,jblk)%icat)
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

            CALL mpi_barrier (p_comm_glb, p_err)

            IF (associated(thisinfo%next)) THEN
               thisinfo => thisinfo%next
            ELSE
               EXIT
            ENDIF

         ENDDO

      ELSEIF (p_is_work) THEN

         CALL mpi_bcast (numblocks, 1, MPI_INTEGER, p_master_address, p_comm_glb, p_err)

         DO WHILE (numblocks > 0)

            mesg(1:2) = (/p_iam_glb, -1/)
            CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_master_address, 0, p_comm_glb, p_err)

            DO WHILE (.true.)

               CALL mpi_recv (catnum, 1, MPI_INTEGER, p_master_address, 1, p_comm_glb, p_stat, p_err)

               IF (catnum == 0) EXIT

               CALL mpi_recv (mesg(1:4), 4, MPI_INTEGER, p_master_address, 1, p_comm_glb, p_stat, p_err)

               imin = mesg(1)
               jmin = mesg(2)
               np   = mesg(3)
               mp   = mesg(4)

               allocate (catch (np,mp))

               CALL mpi_recv (catch, np*mp, MPI_INTEGER,  p_master_address, 1, p_comm_glb, p_stat, p_err)

               nnb = 0
               allocate (nbindex  (np*mp))
               allocate (lenborder(np*mp))
               lenborder(:) = 0

               DO ithis = 1, np
                  DO jthis = 1, mp
                     IF (catch(ithis,jthis) == catnum) THEN

                        DO inext = ithis-1, ithis+1
                           DO jnext = jthis-1, jthis+1
                              IF ((inext >= 1) .and. (inext <= np) .and. (jnext >= 1) .and. (jnext <= mp) &
                                 .and. ((inext /= ithis) .or. (jnext /= jthis))) THEN
                                 IF (catch(inext,jnext) /= catnum) THEN

                                    IF ((catch(inext,jnext) <= 0) .and. (catch(inext,jnext) /= -9)) CYCLE

                                    igthis = ithis + imin - 1
                                    ignext = inext + imin - 1
                                    jgthis = jthis + jmin - 1;  IF (jgthis > mglb) jgthis = jgthis - mglb
                                    jgnext = jnext + jmin - 1;  IF (jgnext > mglb) jgnext = jgnext - mglb

                                    blen = dist_between (igthis, jgthis, ignext, jgnext)

                                    IF (abs(inext-ithis)+abs(jnext-jthis) == 1) THEN
                                       blen = 0.5 * blen
                                    ELSE
                                       blen = 0.25 * blen
                                    ENDIF

                                    IF (nnb > 0) THEN
                                       inb = findloc(nbindex(1:nnb),catch(inext,jnext),dim=1)
                                       IF (inb <= 0) THEN
                                          nnb = nnb + 1
                                          nbindex(nnb) = catch(inext,jnext)
                                          inb = nnb
                                       ENDIF
                                    ELSE
                                       nnb = 1
                                       nbindex(nnb) = catch(inext,jnext)
                                       inb = nnb
                                    ENDIF

                                    lenborder(inb) = lenborder(inb) + blen

                                 ENDIF
                              ENDIF
                           ENDDO
                        ENDDO

                     ENDIF
                  ENDDO
               ENDDO

               IF (nnb == 0) THEN
                  write(*,100) catnum
                  100 format('(S4) Catchment neighbour, ID ', I10, ' has no neighbours.')
               ELSE
                  write(*,101) catnum, nnb, minval(lenborder(1:nnb)), maxval(lenborder(1:nnb))
                  101 format('(S4) Catchment neighbour, ID ', I10, ' has ', I10, ' neighbours (min ', &
                    ES8.2, ', ave ', ES8.2, ')')
               ENDIF

               mesg(1:2) = (/p_iam_glb, catnum/)
               CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_master_address, 0, p_comm_glb, p_err)

               CALL mpi_send (nnb, 1, MPI_INTEGER, p_master_address, 2, p_comm_glb, p_err)
               IF (nnb > 0) THEN
                  CALL mpi_send (nbindex  (1:nnb), nnb, MPI_INTEGER, p_master_address, 2, p_comm_glb, p_err)
                  CALL mpi_send (lenborder(1:nnb), nnb, MPI_REAL4  , p_master_address, 2, p_comm_glb, p_err)
               ENDIF

               deallocate (catch)
               deallocate (nbindex)
               deallocate (lenborder)

            ENDDO

            CALL mpi_barrier (p_comm_glb, p_err)

            numblocks = numblocks - 1

         ENDDO

      ENDIF

   END SUBROUTINE get_basin_neighbour

END MODULE neighbour_mod
