#Makefile include include.mk
############################## Change Log ##################################
# 1.0.0.2
#
# 000908 MJB include.mk-mrc ##
#            Added MAKE environment varaible.
#            Added free format option to F_OPTS for some platforms. ##
# 000907 MJB include.mk-mrc ##
#            Changed the defualts to no NCAR Graphics and no parallel.
#            Also commented out the machine specifics to force the user to
#            select the appropriate machine for them. ##
# 000823 MJB include.mk-mrc ##
#            New - defines all make environment varaibles and is included
#            in all make files. ##
#
############################################################################

# Define make (gnu make works best).

MAKE=/usr/bin/make

# libraries.

BASE=$(ED_ROOT)/build/

# Activate appropriate parts below, comment out others.

#---------------------------------------------------------------
# If you are using a real distribution of NCAR Graphics...
#NCARG_DIR=/usr/local/ncarg-4.3.0/lib
#LIBNCARG=-L$(NCARG_DIR) -lncarg -lncarg_gks -lncarg_c -L/usr/X11R6/lib -lX11 -ldl
#---------------------------------------------------------------

# HDF 5  Libraries
# ED2 HAS OPTIONAL HDF 5 I/O
# If you wish to use this functionality specify USE_HDF5=1
# and specify the location of the include directory
# library files. Make sure you include the zlib.a location too.
USE_HDF5=1
HDF5_DIR=/tigress-hsm/dmedvigy/hdf5-1.8.7/hdf5
HDF5_INCS=-I$(HDF5_DIR)/include -I$(HDF5_DIR)/lib
HDF5_LIBS=-L$(HDF5_DIR)/lib -lhdf5 -lm -lhdf5_fortran -lhdf5 -lhdf5_hl -lz

#---------------------------------------------------------------
# If you have a version of hdf5 compiled in parallel, then you
# may benefit from collective I/O, then use this flag = 1
# Otherwise, set it to zero.

USE_COLLECTIVE_MPIO=0
 
#---------------------------------------------------------------


# netCDF libraries ---------------------------------------------
# If you have netCDF set USENC=1 and type the lib folder
# at NC_LIBS, with the leading -L (e.g. -L/usr/local/lib).
# If you don't have it, leave USENC=0 and type a dummy
# folder for NC_LIBS (e.g. -L/dev/null or leave it blank)
USENC=0
NC_LIBS=-L/dev/null
# --------------------------------------------------------------

# interface ----------------------------------------------------
# This should be 1 unless you are running with -gen-interfaces.
# Interfaces usually make the compilation to crash when the 
# -gen-interfaces option are on, so this flag bypass all 
# interfaces in the code.
USE_INTERF=1


# MPI_Wtime. ---------------------------------------------------
# If USE_MPIWTIME=1, then it will use MPI libraries to compute
# the wall time (the only double-precision intrinsic).  In case
# you don't have it, leave USE_MPIWTIME=0, in which case it will
# use a simpler, single-precision function.
USE_MPIWTIME=1

#----------------- SGI -n32 ------------------------------------
#CMACH=SGI
#F_COMP=f90
#F_OPTS=-n32 -mips4 -r10000 -O2 -freeform \
#       -OPT:IEEE_arithmetic=3:roundoff=3 -OPT:fold_arith_limit=3200
#C_COMP=cc
#C_OPTS=-O -n32 -mips4 -DUNDERSCORE
#LOADER=f90
#LOADER_OPTS=-n32
#LIBS=-lm

#----------------- IBM ------------------------------------
#CMACH=IBM
#F_COMP=f77
#F_OPTS=-O -NA18000 -NQ30000 -qcharlen=10000 -k
#C_COMP=cc
#C_OPTS=-O
#C_PP_OPTS=-P
#LOADER=f77
#LOADER_OPTS=
#LIBS=
#-----------------------------------------------------------------

#----------------- IBM mpxl ------------------------------------
#CMACH=IBM
#F_COMP=mpxlf90
#F_OPTS=-O3 -qstrict -qtune=pwr3 -qarch=pwr3 -NA18000 -NQ30000 -qcharlen=10000 \
#       -qmaxmem=-1
#C_COMP=mpcc
#C_OPTS=-w -O
#C_PP_OPTS=-P
#LOADER=mpxlf90
#LOADER_OPTS=-bmaxdata:1024000000 -bmaxstack:1024000000 -lmassts
#LIBS=
#-----------------------------------------------------------------

#----------------- Sun -------------------------------------------
#CMACH=SUN
#F_COMP=f90
#F_OPTS=-O4
#C_COMP=cc
#C_OPTS=
#LOADER=f90
#LOADER_OPTS=
#LIBS=
#-----------------------------------------------------------------

#----------------- HP/Exemplar ------------------------------------
#CMACH=HP
#F_COMP=f90
#F_OPTS=+O2 +source=free +U77 +Olibcalls +Odataprefetch +Ofltacc
#C_COMP=cc
#C_OPTS=-O
#LOADER=f90
#LOADER_OPTS=
#LIBS=-Wl,-aarchive_shared -lm -lcnx_syscall /lib/libail.sl -lU77 -lc      
#-----------------------------------------------------------------

#----------------- HP/K460 ---------------------------------------
#CMACH=HP
#F_COMP=f90
#F_OPTS=+O2 +U77 +Olibcalls +Odataprefetch +Ofltacc +DA2.0 +DS2.0
#C_COMP=/usr/bin/cc
#C_OPTS=-O -Aa +e -D_HPUX_SOURCE
#LOADER=f90
#LOADER_OPTS=
#LIBS=-lm
#-----------------------------------------------------------------

#----------------------DEC/Compaq Alpha--------------------------
#CMACH=ALPHA
#F_COMP=f90
#F_OPTS=-O -fpe4 -NA18000 -NQ30000 -qcharlen=10000
#C_COMP=cc
#C_OPTS=-DUNDERSCORE -DLITTLE
#LOADER=f90
#LOADER_OPTS=
#LIBS=
#-----------------------------------------------------------------

#-----------------  LINUX Portland Group pgf77/gcc ---------------
#CMACH=PC_LINUX1
#F_COMP=pgf90
#F_OPTS=-Mvect=cachesize:262144,sse -Munroll -Mnoframe -O2 -pc 64 \
#       -Mfree
#  -Mbyteswapio

# -tp athlonxp -fastsse

#C_COMP=gcc
#C_COMP=pgcc
#C_OPTS=-O3 -DUNDERSCORE -DLITTLE
#LOADER=pgf90
#LOADER_OPTS=-v -Wl,
#LIBS=
#-----------------------------------------------------------------

#----------------  NEC-SX6 ---------------
#CMACH=NEC_SX
#F_COMP=sxmpif90
#C_COMP=sxmpic++
#LOADER=sxmpif90
#C_LOADER=sxmpic++
#LIBS=
#MOD_EXT=mod
#Compiler options:
#F_OPTS=-eC -C debug -Wf "-cont -init stack=nan heap=nan"
#F_OPTS=-ftrace
#F_OPTS=-C ssafe
#F_OPTS=
#C_OPTS=-f90lib
#C_OPTS=-f90lib -C ssafe
#C_OPTS=-f90lib -ftrace
#LOADER_OPTS=-eC -C debug -f90lib -Wf "-cont -init stack=nan heap=nan"
#LOADER_OPTS=-C debug -f90lib
#LOADER_OPTS=-ftrace
#LOADER_OPTS=
#-----------------------------------------------------------------

#-----------------  LINUX INTEL FORTRAN-95 Compiler/GCC  ---------
CMACH=PC_LINUX1
F_COMP=mpif90
C_COMP=gcc
LOADER=mpif90
C_LOADER=mpicc
LIBS=
MOD_EXT=mod

##################################### COMPILER OPTIONS #####################################
#------------------------------------------------------------------------------------------#
# A. Pickiest - Use this whenever you change arguments on functions and subroutines.       #
#               This will perform the same tests as B but it will also check whether all   #
#               arguments match between subroutine declaration and subroutine calls.       #
#               WARNING: In order to really check all interfaces you must compile with     #
#                        this option twice:                                                #
#               1. Compile (./compile.sh)                                                  #
#               2. Prepare second compilation(./2ndcomp.sh)                                #
#               3. Compile one more time (./compile.sh)                                    #
#               If the compilation fails either at step 1 or 3, then your code has inter-  #
#                  face problems. If it successfully compiles, then you can switch to B.   #
# B. Pickiest with no interface - This will compile fast but the run will be slow due to   #
#    the -O0 option. However, by setting -O0 you will take full advantage of the intel     #
#    debugger.                                                                             #
#    Ideally, you should compile your code with this option whenever you make any changes. #
#    Note, however, that if you change arguments you should first try A.                   #
# C. Fast debugging - This will check pretty much the same as B, but with higher optimiza- #
#    tion. However, by setting -O2 you won't be able to see all variables on the debugger. #
#    You should use this if your run seems to be working fine at the beginning, but it is  #
#    failing or giving instabilities, or funny results after a long time.                  #
# D. Fast check - This will check pretty much the same as C, but it will not set up        #
#    anything for the debugger. Use this only if you really don't want to deal with idb or #
#    if you have a good idea of which problem you are dealing with.                        #
# E. Fast - This is all about performance, use only when you are sure that the model has   #
#           no code problem, and you want results asap. This will not check for any        #
#           problems, which means that this is an option suitable for end users, not de-   #
#           velopers.                                                                      #
#------------------------------------------------------------------------------------------#
KIND_COMP=E
#E

#------------------------------------------------------------------------------------------#
ifeq ($(KIND_COMP),A)
   USE_INTERF=0
   F_OPTS= -FR -O0 -recursive -Vaxlib  -check all -g -fpe0 -ftz -gen-interfaces            \
            -warn interfaces -debug extended -debug inline_debug_info                      \
             -debug-parameters all -traceback -ftrapuv
   C_OPTS= -O0 -DLITTLE  -g -traceback -debug extended
   LOADER_OPTS= -FR -O0 -Vaxlib  -check all -g -fpe0 -ftz -gen-interfaces \
              -warn interfaces -debug extended -debug inline_debug_info  \
              -debug-parameters all -traceback -ftrapuv
   C_LOADER_OPTS=-v -g -traceback
   #---------------------------------------------------------------------------------------#
endif
ifeq ($(KIND_COMP),B)
   USE_INTERF=1
   F_OPTS= -FR -O0 -recursive -Vaxlib  -check all -g -fpe0 -ftz  -debug extended           \
           -debug inline_debug_info -debug-parameters all -traceback -ftrapuv
   C_OPTS= -O0 -DLITTLE  -g -traceback -debug extended
   LOADER_OPTS= -FR -O0 -Vaxlib  -check all -g -fpe0 -ftz -debug extended                  \
                -debug inline_debug_info -debug-parameters all -traceback -ftrapuv
   C_LOADER_OPTS=-v -g -traceback
   #---------------------------------------------------------------------------------------#
endif
ifeq ($(KIND_COMP),C)
    USE_INTERF=1
    F_OPTS= -FR -O2 -recursive -Vaxlib  -check all -g -fpe0 -ftz  -debug extended          \
            -debug inline_debug_info -debug-parameters all -traceback -ftrapuv
    C_OPTS= -O2 -DLITTLE  -g -traceback -debug extended
    LOADER_OPTS= -FR -O2 -Vaxlib  -check all -g -fpe0 -ftz -debug extended                 \
                 -debug inline_debug_info -debug-parameters all -traceback -ftrapuv
    C_LOADER_OPTS=-v -g -traceback
   #---------------------------------------------------------------------------------------#
endif
ifeq ($(KIND_COMP),D)
   USE_INTERF=1
   F_OPTS= -FR -O2 -recursive -Vaxlib -check all -fpe0 -ftz  -traceback -ftrapuv
   C_OPTS= -O2 -DLITTLE -traceback
   LOADER_OPTS= -FR -O2 -Vaxlib  -check all -fpe0 -ftz -traceback -ftrapuv
   C_LOADER_OPTS=-v -traceback
   #---------------------------------------------------------------------------------------#
endif
ifeq ($(KIND_COMP),E)
   USE_INTERF=1
   F_OPTS=-xT -O3 -fno-alias -free -traceback
#   F_OPTS=-g -fp-model precise -check bounds -traceback \
#          -debug extended -check uninit -ftrapuv
   C_OPTS= -O3 -DUNDERSCORE -DLITTLE
   LOADER_OPTS=-i-static $(F_OPTS)
   C_LOADER_OPTS=-v -traceback


   #---------------------------------------------------------------------------------------#
endif
#------------------------------------------------------------------------------------------#

############################################################################################

#-----------------------------------------------------------------

# If compiling for a single-CPU platform only (without MPI):

#-----------------------------------------------------------------
#PAR_LIBS=
#PAR_DEFS=
#-----------------------------------------------------------------

# Else if using MPI libraries:

#---------------SGI-----------------------------------------------
#with mpich parallel stuff
#MPI_PATH=/n/Moorcroft_Lab/Users/mlongo/util/mpich
#PAR_INCS=-I$(MPI_PATH)/include
#PAR_LIBS=-L$(MPI_PATH)/lib/IRIXN32/ch_shmem -lmpi
#  or with SGI Parallel stuff
#PAR_LIBS=-L/usr/lib32 -lmpi
#  need this for both
#PAR_DEFS=-DRAMS_MPI
#-----------------------------------------------------------------

#---------------IBM-----------------------------------------------
#MPI_PATH=/usr/local/mpich
#PAR_INCS=-I$(MPI_PATH)/include
#PAR_LIBS=-L$(MPI_PATH)/lib/rs6000/ch_p4 -lmpi 
#PAR_DEFS=-DRAMS_MPI
#-----------------------------------------------------------------

#---------------Sun-----------------------------------------------
#MPI_PATH=/usr/local/mpich
#PAR_INCS=-I$(MPI_PATH)/include
#PAR_LIBS=-L$(MPI_PATH)/lib/solaris/ch_p4 -lmpi 
#PAR_DEFS=-DRAMS_MPI
#-----------------------------------------------------------------

#---------------HP-Exemplar---------------------------------------
#MPI_PATH=/opt/mpi
#PAR_INCS=-I$(MPI_PATH)/include
#PAR_LIBS=$(MPI_PATH)/lib/pa1.1/libmpi.a
#PAR_DEFS=-DRAMS_MPI
#-----------------------------------------------------------------

#---------------LINUX Portland Group pgf77/gcc--------------------
#MPI_PATH=/n/Moorcroft_Lab/Lab/apps/i91
#PAR_INCS=-I$(MPI_PATH)/include
#PAR_LIBS=-L$(MPI_PATH)/lib -lmpich -lpmpich
#PAR_DEFS=-DRAMS_MPI
#-----------------------------------------------------------------

#---------------If using scritps 'mpicc' e 'mpif90'---------------'
MPI_PATH=
PAR_INCS=
PAR_LIBS=
PAR_DEFS=-DRAMS_MPI
#-----------------------------------------------------------------

# For IBM,HP,SGI,ALPHA,LINUX use these:
ARCHIVE=ar rs
# For NEC SX-6
#ARCHIVE=sxar rs
# For SUN,CONVEX
#ARCHIVE=ar r'

