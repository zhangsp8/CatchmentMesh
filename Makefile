
FF = /usr/bin/mpif90
FOPTS = -fdefault-real-8 -ffree-form -C -g -u -xcheck=stkovf \
		-ffpe-trap=invalid,zero,overflow,underflow -fbounds-check \
		-mcmodel=medium -fbacktrace -fdump-core -cpp -fopenmp
NETCDF_LIB = /usr/lib/x86_64-linux-gnu
NETCDF_INC = /usr/include

# FF = mpif90
# FOPTS = -r8 -free -g -check uninit -check bounds -check pointers -traceback
# NETCDF_INC = /opt/netcdf/include
# NETCDF_LIB = /opt/netcdf/lib

objs = precision.o     \
		 MathConstants.o \
		 utils_mod.o     \
		 cvt_mod.o \
		 task_mod.o  \
		 ncio_serial.o \
		 hydro_data_mod.o \
		 river_lake_mod.o \
		 catchment_mod.o \
		 hillslope_mod.o \
		 information_mod.o \
		 output_mod.o

main : $(objs) main.F90
	$(FF) -c $(FOPTS) -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdff -o $@.o $@.F90
	$(FF) $(FOPTS) $(objs) $@.o -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdff -o $@ 

get_rivermouth : $(objs) get_rivermouth.F90
	$(FF) -c $(FOPTS) -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdff -o $@.o $@.F90
	$(FF) $(FOPTS) $(objs) $@.o -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdff -o $@ 

$(objs) : %.o : %.F90
	$(FF) -c $(FOPTS) -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdff -o $@ $<

.PHONY : clean

clean : 
	-rm -f *.o *.mod
