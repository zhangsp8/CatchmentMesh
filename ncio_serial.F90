MODULE ncio_serial

   USE netcdf
   IMPLICIT NONE

   ! PUBLIC subroutines

   PUBLIC :: ncio_create_file
   PUBLIC :: ncio_inquire_length
   PUBLIC :: ncio_define_dimension

   interface ncio_read_serial
      MODULE procedure ncio_read_serial_real8_1d 
      MODULE procedure ncio_read_serial_int8_2d 
      MODULE procedure ncio_read_serial_int32_2d 
      MODULE procedure ncio_read_serial_real4_2d 
   END interface ncio_read_serial

   interface ncio_write_serial
      MODULE procedure ncio_write_serial_int32_0d 
      MODULE procedure ncio_write_serial_int32_1d 
      MODULE procedure ncio_write_serial_real8_1d 
      MODULE procedure ncio_write_serial_real4_1d 
      MODULE procedure ncio_write_serial_int8_2d 
      MODULE procedure ncio_write_serial_int32_2d 
      MODULE procedure ncio_write_serial_real4_2d 
   END interface ncio_write_serial

CONTAINS

   ! ----
   SUBROUTINE nccheck (status)
      integer, intent(in) :: status

      IF (status /= NF90_NOERR) THEN
         print *, trim(nf90_strerror(status))
         STOP 2
      ENDIF
   END SUBROUTINE nccheck

   ! ----
   SUBROUTINE ncio_create_file (filename)

      USE netcdf
      IMPLICIT NONE

      character(len=*), intent(in)  :: filename

      ! Local Variables
      integer :: ncid

      CALL nccheck( nf90_create(trim(filename), ior(NF90_CLOBBER,NF90_NETCDF4), ncid) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_create_file

   !---------------------------------------------------------
   SUBROUTINE ncio_inquire_varsize (filename, dataname, varsize)

      USE netcdf
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      integer, allocatable, intent(out) :: varsize(:)

      ! Local variables
      integer :: ncid, varid, ndims, idm
      integer, allocatable :: dimids(:)

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )

      CALL nccheck( nf90_inquire_variable(ncid, varid, ndims = ndims) )
      allocate (dimids(ndims))
      CALL nccheck( nf90_inquire_variable(ncid, varid, dimids = dimids) )

      allocate (varsize(ndims))
      DO idm = 1, ndims
         CALL nccheck( nf90_inquire_dimension(ncid, dimids(idm), len = varsize(idm)) )
      ENDDO 
      
      CALL nccheck( nf90_close(ncid) )
      deallocate (dimids)

   END SUBROUTINE ncio_inquire_varsize

   !---------------------------------------------------------
   SUBROUTINE ncio_inquire_length (filename, dataname, length)

      USE netcdf
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      integer, intent(out) :: length

      ! Local variables
      integer :: ncid, varid, ndims
      integer, allocatable :: dimids(:)

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )

      CALL nccheck( nf90_inquire_variable(ncid, varid, ndims = ndims) )
      allocate (dimids(ndims))
      CALL nccheck( nf90_inquire_variable(ncid, varid, dimids = dimids) )
      CALL nccheck( nf90_inquire_dimension(ncid, dimids(ndims), len = length) )
      
      CALL nccheck( nf90_close(ncid) )
      deallocate (dimids)

   END SUBROUTINE ncio_inquire_length

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_real8_1d (filename, dataname, rdata)

      USE netcdf
      USE precision
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      real(r8), allocatable, intent(out) :: rdata (:)

      ! Local variables
      integer :: ncid, varid
      integer, allocatable :: varsize(:)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_real8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_int8_2d (filename, dataname, rdata)

      USE netcdf
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      integer(1), allocatable, intent(out) :: rdata (:,:)

      ! Local variables
      integer :: ncid, varid
      integer, allocatable :: varsize(:)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_int8_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_int32_2d (filename, dataname, rdata)

      USE netcdf
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      integer(4), allocatable, intent(out) :: rdata (:,:)

      ! Local variables
      integer :: ncid, varid
      integer, allocatable :: varsize(:)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_int32_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_real4_2d (filename, dataname, rdata)

      USE netcdf
      USE precision
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      real(4), allocatable, intent(out) :: rdata (:,:)

      ! Local variables
      integer :: ncid, varid
      integer, allocatable :: varsize(:)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_real4_2d

   ! -------------------------------
   SUBROUTINE ncio_define_dimension (filename, dimname, dimlen)

      USE netcdf
      USE precision
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dimname
      integer, intent(in) :: dimlen

      ! Local variables
      integer :: ncid, dimid, status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )

      status = nf90_inq_dimid(ncid, trim(dimname), dimid)
      IF (status /= NF90_NOERR) THEN
         CALL nccheck (nf90_redef(ncid))
         IF (dimlen == 0) THEN
            CALL nccheck( nf90_def_dim(ncid, trim(dimname), NF90_UNLIMITED, dimid) )
         ELSE
            CALL nccheck( nf90_def_dim(ncid, trim(dimname), dimlen, dimid) )
         ENDIF 
         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_define_dimension

   ! ----
   SUBROUTINE ncio_put_attr_str (filename, varname, attrname, attrval)

   USE netcdf
   IMPLICIT NONE

   character(len=*), intent(in)  :: filename, varname, attrname, attrval

   ! Local Variables
   integer :: ncid, varid

      CALL nccheck( nf90_open (trim(filename), NF90_WRITE, ncid) )
      CALL nccheck (nf90_inq_varid (ncid, trim(varname), varid))
      CALL nccheck (nf90_redef (ncid))
      CALL nccheck (nf90_put_att (ncid, varid, trim(attrname), trim(attrval)))
      CALL nccheck (nf90_enddef (ncid))
      CALL nccheck( nf90_close (ncid))

   END SUBROUTINE ncio_put_attr_str

   ! ----
   SUBROUTINE ncio_put_attr_real4 (filename, varname, attrname, attrval)

   USE netcdf
   IMPLICIT NONE

   character(len=*), intent(in)  :: filename, varname, attrname
   real(kind = 4),   intent(in)  :: attrval

   ! Local Variables
   integer :: ncid, varid

      CALL nccheck( nf90_open (trim(filename), NF90_WRITE, ncid) )
      CALL nccheck (nf90_inq_varid (ncid, trim(varname), varid))
      CALL nccheck (nf90_redef (ncid))
      CALL nccheck (nf90_put_att (ncid, varid, trim(attrname), attrval))
      CALL nccheck (nf90_enddef (ncid))
      CALL nccheck( nf90_close (ncid))

   END SUBROUTINE ncio_put_attr_real4

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int32_0d (filename, dataname, wdata)

      USE netcdf
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: wdata 

      ! Local variables
      integer :: ncid, varid, status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, varid = varid))
         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int32_0d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int32_1d (filename, dataname, wdata, dimname, compress)

      USE netcdf
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: wdata (:)

      character(len=*), intent(in), optional :: dimname
      integer, intent(in), optional :: compress

      ! Local variables
      integer :: ncid, varid, dimid, status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. present(dimname)) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dimname), dimid))
         
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int32_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_1d (filename, dataname, wdata, dimname, compress)

      USE netcdf
      USE precision
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      real(r8), intent(in) :: wdata (:)

      character(len=*), intent(in), optional :: dimname
      integer, intent(in), optional :: compress

      ! Local variables
      integer :: ncid, varid, dimid, status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. present(dimname)) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dimname), dimid))
         
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real4_1d (filename, dataname, wdata, dimname, compress)

      USE netcdf
      USE precision
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      real(kind = 4), intent(in) :: wdata (:)

      character(len=*), intent(in), optional :: dimname
      integer, intent(in), optional :: compress

      ! Local variables
      integer :: ncid, varid, dimid, status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. present(dimname)) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dimname), dimid))
         
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_FLOAT, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_FLOAT, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real4_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int8_2d (filename, dataname, wdata, &
         dim1name, dim2name, compress)

      USE netcdf
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      integer(1), intent(in) :: wdata (:,:)
      
      character(len=*), intent(in), optional :: dim1name, dim2name
      integer, intent(in), optional :: compress

      ! Local variables
      integer :: ncid, varid, dimid(2), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_BYTE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_BYTE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int8_2d
   
   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int16_2d (filename, dataname, wdata, &
         dim1name, dim2name, compress)

      USE netcdf
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      integer(2), intent(in) :: wdata (:,:)
      
      character(len=*), intent(in), optional :: dim1name, dim2name
      integer, intent(in), optional :: compress

      ! Local variables
      integer :: ncid, varid, dimid(2), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_SHORT, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_SHORT, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int16_2d
   
   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int32_2d (filename, dataname, wdata, &
         dim1name, dim2name, compress)

      USE netcdf
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: wdata (:,:)
      
      character(len=*), intent(in), optional :: dim1name, dim2name
      integer, intent(in), optional :: compress

      ! Local variables
      integer :: ncid, varid, dimid(2), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int32_2d
   
   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real4_2d (filename, dataname, wdata, &
         dim1name, dim2name, compress)

      USE netcdf
      USE precision
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      real(4), intent(in) :: wdata (:,:)
      
      character(len=*), intent(in), optional :: dim1name, dim2name
      integer, intent(in), optional :: compress

      ! Local variables
      integer :: ncid, varid, dimid(2), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_FLOAT, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_FLOAT, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real4_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_2d (filename, dataname, wdata, &
         dim1name, dim2name, compress)

      USE netcdf
      USE precision
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      real(r8), intent(in) :: wdata (:,:)
      
      character(len=*), intent(in), optional :: dim1name, dim2name
      integer, intent(in), optional :: compress

      ! Local variables
      integer :: ncid, varid, dimid(2), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int32_3d (filename, dataname, wdata, &
         dim1name, dim2name, dim3name, compress)

      USE netcdf
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      integer, intent(in) :: wdata (:,:,:)
      
      character(len=*), intent(in), optional :: dim1name, dim2name, dim3name
      integer, intent(in), optional :: compress

      ! Local variables
      integer :: ncid, varid, dimid(3), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name) .and. present(dim3name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim3name), dimid(3)))
         
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int32_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_3d (filename, dataname, wdata, &
         dim1name, dim2name, dim3name, compress)

      USE netcdf
      USE precision
      IMPLICIT NONE

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: dataname
      real(r8), intent(in) :: wdata (:,:,:)
      
      character(len=*), intent(in), optional :: dim1name, dim2name, dim3name
      integer, intent(in), optional :: compress

      ! Local variables
      integer :: ncid, varid, dimid(3), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name) .and. present(dim3name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim3name), dimid(3)))
         
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_3d


END MODULE ncio_serial
