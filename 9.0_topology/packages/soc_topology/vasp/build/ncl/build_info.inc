    character(len=*), parameter :: cpp_options = '&
&-DHOST="LinuxIFC" &
&-DMPI &
&-DMPI_BLOCK=8000 &
&-Duse_collective &
&-DscaLAPACK &
&-DCACHE_SIZE=4000 &
&-Davoidalloc &
&-Dvasp6 &
&-Dtbdyn &
&-Dfock_dblbuf &
&-DVASP_HDF5 &
&-DVASP2WANNIER90'
    character(len=*), parameter :: link_line   = '&
&mpiifort &
&-fc=ifx &
&-mkl=sequential &
&-Llib &
&-ldmy &
&-Lparser &
&-lparser &
&-lstdc++ &
&-L/usr/physics/oneapi/2024.2.1/mkl/2024.2/lib/intel64 &
&-lmkl_scalapack_lp64 &
&-lmkl_blacs_intelmpi_lp64 &
&-L/usr/physics/hdf/5/1.14.1-2_intel2021/lib &
&-lhdf5_fortran &
&/cmt2/lniu6305/Packages/soc_topology_20260922_1435/wannier90-3.1.0/libwannier.a'
    character(len=*), parameter :: fc     = '&
&mpiifort &
&-fc=ifx'
    character(len=*), parameter :: fcl    = '&
&mpiifort &
&-fc=ifx &
&-mkl=sequential'
    character(len=*), parameter :: fflags = '&
&-assume &
&byterecl &
&-w &
&-axCORE-AVX512'
    character(len=*), parameter :: llibs  = '&
&-lstdc++ &
&-L/usr/physics/oneapi/2024.2.1/mkl/2024.2/lib/intel64 &
&-lmkl_scalapack_lp64 &
&-lmkl_blacs_intelmpi_lp64 &
&-L/usr/physics/hdf/5/1.14.1-2_intel2021/lib &
&-lhdf5_fortran &
&/cmt2/lniu6305/Packages/soc_topology_20260922_1435/wannier90-3.1.0/libwannier.a'
    character(len=*), parameter :: incs   = '&
&-I/usr/physics/oneapi/2024.2.1/mkl/2024.2/include/fftw &
&-I/usr/physics/hdf/5/1.14.1-2_intel2021/include'
