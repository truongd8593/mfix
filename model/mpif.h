! "USMID @(#) mpi_t3e/MPI.inc	11.7	01/05/2000 18:49:31"
!
!	(C) COPYRIGHT SILICON GRAPHICS, INC.
!	UNPUBLISHED PROPRIETARY INFORMATION.
!	ALL RIGHTS RESERVED.
!
!
! ------------------------------
! This is mpif.h for T3E systems
! ------------------------------
!
! Copyright Notice
!
! (c) Copyright 1995, 1996, 1997, The University of Edinburgh.
!
! ======================================================================
!									
!     File: MPI.inc							
!  Project: KTG-T3DMPI							
!   Author: K.Cameron & A.G. Smith & K.J. Wierenga			
!  Created: 07/11/1994							
!  Descrip: The MPI Fortran interface header file			
! 									
! ======================================================================

! === Follows MPI (Message-passing Interface) Standard   Annex A. ===

! Version of MPI standard: VERSION.SUBVERSION 
       integer    MPI_VERSION, MPI_SUBVERSION
       parameter (MPI_VERSION    = 1)
       parameter (MPI_SUBVERSION = 2)

! MPI-1 Return codes 
       integer MPI_SUCCESS
       parameter (MPI_SUCCESS			= 0)

       integer MPI_ERR_BUFFER
       integer MPI_ERR_COUNT
       integer MPI_ERR_TYPE
       integer MPI_ERR_TAG
       integer MPI_ERR_COMM
       integer MPI_ERR_RANK
       integer MPI_ERR_REQUEST
       integer MPI_ERR_ROOT
       integer MPI_ERR_GROUP
       integer MPI_ERR_OP
       integer MPI_ERR_TOPOLOGY
       integer MPI_ERR_DIMS
       integer MPI_ERR_ARG
       integer MPI_ERR_UNKNOWN
       integer MPI_ERR_TRUNCATE
       integer MPI_ERR_OTHER
       integer MPI_ERR_INTERN
       integer MPI_ERR_IN_STATUS
       integer MPI_ERR_PENDING

       parameter (MPI_ERR_BUFFER                = 1)
       parameter (MPI_ERR_COUNT                 = 2)
       parameter (MPI_ERR_TYPE                  = 3)
       parameter (MPI_ERR_TAG                   = 4)
       parameter (MPI_ERR_COMM                  = 5)
       parameter (MPI_ERR_RANK                  = 6)
       parameter (MPI_ERR_REQUEST               = 7)
       parameter (MPI_ERR_ROOT                  = 8)
       parameter (MPI_ERR_GROUP                 = 9)
       parameter (MPI_ERR_OP                    = 10)
       parameter (MPI_ERR_TOPOLOGY              = 11)
       parameter (MPI_ERR_DIMS                  = 12)
       parameter (MPI_ERR_ARG                   = 13)
       parameter (MPI_ERR_UNKNOWN               = 14)
       parameter (MPI_ERR_TRUNCATE              = 15)
       parameter (MPI_ERR_OTHER                 = 16)
       parameter (MPI_ERR_INTERN                = 17)
       parameter (MPI_ERR_IN_STATUS             = 18)
       parameter (MPI_ERR_PENDING               = 19)

! MPI-2 Return codes

       integer MPI_ERR_ACCESS
       integer MPI_ERR_AMODE
       integer MPI_ERR_ASSERT
       integer MPI_ERR_BAD_FILE
       integer MPI_ERR_BASE
       integer MPI_ERR_CONVERSION
       integer MPI_ERR_DISP
       integer MPI_ERR_DUP_DATAREP
       integer MPI_ERR_FILE_EXISTS
       integer MPI_ERR_FILE_IN_USE
       integer MPI_ERR_FILE
       integer MPI_ERR_INFO_KEY
       integer MPI_ERR_INFO_NOKEY
       integer MPI_ERR_INFO_VALUE
       integer MPI_ERR_INFO
       integer MPI_ERR_IO
       integer MPI_ERR_KEYVAL
       integer MPI_ERR_LOCKTYPE
       integer MPI_ERR_NAME
       integer MPI_ERR_NO_MEM
       integer MPI_ERR_NOT_SAME
       integer MPI_ERR_NO_SPACE
       integer MPI_ERR_NO_SUCH_FILE
       integer MPI_ERR_PORT
       integer MPI_ERR_QUOTA
       integer MPI_ERR_READ_ONLY
       integer MPI_ERR_RMA_CONFLICT
       integer MPI_ERR_RMA_SYNC
       integer MPI_ERR_SERVICE
       integer MPI_ERR_SIZE
       integer MPI_ERR_SPAWN
       integer MPI_ERR_UNSUPPORTED_DATAREP
       integer MPI_ERR_UNSUPPORTED_OPERATION
       integer MPI_ERR_WIN

       parameter (MPI_ERR_ACCESS                  = 28)
       parameter (MPI_ERR_AMODE                   = 29)
       parameter (MPI_ERR_ASSERT                  = 30)
       parameter (MPI_ERR_BAD_FILE                = 31)
       parameter (MPI_ERR_BASE                    = 32)
       parameter (MPI_ERR_CONVERSION              = 33)
       parameter (MPI_ERR_DISP                    = 34)
       parameter (MPI_ERR_DUP_DATAREP             = 35)
       parameter (MPI_ERR_FILE_EXISTS             = 36)
       parameter (MPI_ERR_FILE_IN_USE             = 37)
       parameter (MPI_ERR_FILE                    = 38)
       parameter (MPI_ERR_INFO_KEY                = 39)
       parameter (MPI_ERR_INFO_NOKEY              = 40)
       parameter (MPI_ERR_INFO_VALUE              = 41)
       parameter (MPI_ERR_INFO                    = 42)
       parameter (MPI_ERR_IO                      = 43)
       parameter (MPI_ERR_KEYVAL                  = 44)
       parameter (MPI_ERR_LOCKTYPE                = 45)
       parameter (MPI_ERR_NAME                    = 46)
       parameter (MPI_ERR_NO_MEM                  = 47)
       parameter (MPI_ERR_NOT_SAME                = 48)
       parameter (MPI_ERR_NO_SPACE                = 49)
       parameter (MPI_ERR_NO_SUCH_FILE            = 50)
       parameter (MPI_ERR_PORT                    = 51)
       parameter (MPI_ERR_QUOTA                   = 52)
       parameter (MPI_ERR_READ_ONLY               = 53)
       parameter (MPI_ERR_RMA_CONFLICT            = 54)
       parameter (MPI_ERR_RMA_SYNC                = 55)
       parameter (MPI_ERR_SERVICE                 = 56)
       parameter (MPI_ERR_SIZE                    = 57)
       parameter (MPI_ERR_SPAWN                   = 58)
       parameter (MPI_ERR_UNSUPPORTED_DATAREP     = 59)
       parameter (MPI_ERR_UNSUPPORTED_OPERATION   = 60)
       parameter (MPI_ERR_WIN                     = 61)

       integer MPI_ERR_LASTCODE
       parameter (MPI_ERR_LASTCODE                = 100)

! Assorted constants 
!     [ MPI_BOTTOM must be recognised by its address:
!       it must reside in a common block, value is not used ]
      integer MPI_BOTTOM
      common /MPIPRIVBOT/ MPI_BOTTOM
      integer    MPI_PROC_NULL
      parameter (MPI_PROC_NULL       = -1)
      integer    MPI_ANY_SOURCE
      parameter (MPI_ANY_SOURCE      = -101)
      integer    MPI_ANY_TAG
      parameter (MPI_ANY_TAG         = -102)
      integer    MPI_UNDEFINED
      parameter (MPI_UNDEFINED       = -1)
!     MPI_UB -- as elementary datatype below
!     MPI_LB -- as elementary datatype below
      integer    MPI_BSEND_OVERHEAD
      parameter (MPI_BSEND_OVERHEAD  = 0)
      integer    MPI_KEYVAL_INVALID
      parameter (MPI_KEYVAL_INVALID  = -1)

! --- Status size and reserved index values (Fortran)
      integer    MPI_STATUS_SIZE
      parameter (MPI_STATUS_SIZE     = 6)
      integer    MPI_SOURCE
      parameter (MPI_SOURCE          = 1)
      integer    MPI_TAG
      parameter (MPI_TAG             = 2)
      integer    MPI_ERROR
      parameter (MPI_ERROR           = 3)

! --- MPI-2 STATUS_IGNORE
!     [ MPI_STATUS(ES)_IGNORE must be recognised by their addresses:
!       they must reside in a common block, values are not used.
!     ]
      integer    MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
      integer    MPI_STATUSES_IGNORE(MPI_STATUS_SIZE, 1)
      common /MPIPRIVSTAT/ MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE

! --- Error-handling specifiers
      integer    MPI_ERRORS_ARE_FATAL
      parameter (MPI_ERRORS_ARE_FATAL = 0)
      integer    MPI_ERRORS_RETURN
      parameter (MPI_ERRORS_RETURN    = 1)

! --- Maximum sizes for strings
      integer    MPI_MAX_PROCESSOR_NAME
      parameter (MPI_MAX_PROCESSOR_NAME = 256)
      integer    MPI_MAX_ERROR_STRING
      parameter (MPI_MAX_ERROR_STRING   = 256)

! --- Elementary datatypes (Fortran)
      integer    MPI_INTEGER
      parameter (MPI_INTEGER            = 14)
      integer    MPI_REAL
      parameter (MPI_REAL               = 15)
      integer    MPI_DOUBLE_PRECISION
      parameter (MPI_DOUBLE_PRECISION   = 16)
      integer    MPI_COMPLEX
      parameter (MPI_COMPLEX            = 17)
! ---   [ nb: DOUBLE COMPLEX is optional basic datatype (v1.2) ]
      integer    MPI_LOGICAL
      parameter (MPI_LOGICAL            = 19)
      integer    MPI_CHARACTER
      parameter (MPI_CHARACTER          = 20)
      integer    MPI_BYTE
      parameter (MPI_BYTE               = 21)
      integer    MPI_PACKED
      parameter (MPI_PACKED             = 22)
      integer    MPI_UB
      parameter (MPI_UB                 = 29)
      integer    MPI_LB
      parameter (MPI_LB                 = 30)

! --- Optional datatypes (Fortran)
      integer    MPI_DOUBLE_COMPLEX
      parameter (MPI_DOUBLE_COMPLEX     = 18)
      integer    MPI_INTEGER1
      parameter (MPI_INTEGER1           = 23)
      integer    MPI_INTEGER2
      parameter (MPI_INTEGER2           = 24)
      integer    MPI_INTEGER4
      parameter (MPI_INTEGER4           = 25)
      integer    MPI_INTEGER8
      parameter (MPI_INTEGER8           = 14)
      integer    MPI_REAL4
      parameter (MPI_REAL4              = 27)
      integer    MPI_REAL8
      parameter (MPI_REAL8              = 28)

      integer    MPI_LOGICAL1
      parameter (MPI_LOGICAL1           = 40)
      integer    MPI_LOGICAL2
      parameter (MPI_LOGICAL2           = 41)
      integer    MPI_LOGICAL4
      parameter (MPI_LOGICAL4           = 42)
      integer    MPI_COMPLEX8
      parameter (MPI_COMPLEX8           = 43)
      integer    MPI_COMPLEX16
      parameter (MPI_COMPLEX16          = 44)

! --- Datatypes for reduction functions (Fortran)

      integer    MPI_2REAL
      parameter (MPI_2REAL               = 37)
      integer    MPI_2DOUBLE_PRECISION
      parameter (MPI_2DOUBLE_PRECISION   = 38)
      integer    MPI_2INTEGER
      parameter (MPI_2INTEGER            = 39)
!       MPI_2COMPLEX not defined [MPI Errata Oct 94]

! Reserved communicators 
       integer    MPI_COMM_WORLD
       parameter (MPI_COMM_WORLD      = 0)
       integer    MPI_COMM_SELF
       parameter (MPI_COMM_SELF       = 1)

! Results of communicator and group comparisons 
       integer    MPI_UNEQUAL
       parameter (MPI_UNEQUAL         = -1)
       integer    MPI_SIMILAR
       parameter (MPI_SIMILAR         = -2)
       integer    MPI_IDENT
       parameter (MPI_IDENT           = -3)
       integer    MPI_CONGRUENT
       parameter (MPI_CONGRUENT       = -4)

! Environmental inquiry keys 
       integer    MPI_TAG_UB
       parameter (MPI_TAG_UB          = -50)
       integer    MPI_IO
       parameter (MPI_IO              = -49)
       integer    MPI_HOST
       parameter (MPI_HOST            = -48)
       integer    MPI_WTIME_IS_GLOBAL
       parameter (MPI_WTIME_IS_GLOBAL = -47)

! Collective operations 
       integer    MPI_MAX
       parameter (MPI_MAX             = 0)
       integer    MPI_MIN
       parameter (MPI_MIN             = 1)
       integer    MPI_SUM
       parameter (MPI_SUM             = 2)
       integer    MPI_PROD
       parameter (MPI_PROD            = 3)
       integer    MPI_MAXLOC
       parameter (MPI_MAXLOC          = 4)
       integer    MPI_MINLOC
       parameter (MPI_MINLOC          = 5)
       integer    MPI_BAND
       parameter (MPI_BAND            = 6)
       integer    MPI_BOR
       parameter (MPI_BOR             = 7)
       integer    MPI_BXOR
       parameter (MPI_BXOR            = 8)
       integer    MPI_LAND
       parameter (MPI_LAND            = 9)
       integer    MPI_LOR
       parameter (MPI_LOR             = 10)
       integer    MPI_LXOR
       parameter (MPI_LXOR            = 11)

! MPI-2 collective operations
       integer    MPI_REPLACE
       parameter (MPI_REPLACE         = 12)

! Null handles 
       integer    MPI_GROUP_NULL
       parameter (MPI_GROUP_NULL      = -1)
       integer    MPI_COMM_NULL
       parameter (MPI_COMM_NULL       = -1)
       integer    MPI_DATATYPE_NULL
       parameter (MPI_DATATYPE_NULL   = -1)
       integer    MPI_REQUEST_NULL
       parameter (MPI_REQUEST_NULL    = 0)
       integer    MPI_OP_NULL
       parameter (MPI_OP_NULL         = -1)
       integer    MPI_ERRHANDLER_NULL
       parameter (MPI_ERRHANDLER_NULL = -1)

! MPI-2 null handles
       integer    MPI_WIN_NULL
       parameter (MPI_WIN_NULL        = 0)

!     [ Predefined attribute callbacks must be recognised by address:
!       they must reside in a common block, values are not used ]
      integer    MPI_NULL_COPY_FN
      integer    MPI_DUP_FN
      integer    MPI_NULL_DELETE_FN
      common /MPIPRIVATTRCOPY/ MPI_NULL_COPY_FN,  MPI_DUP_FN
      common /MPIPRIVATTRDEL/  MPI_NULL_DELETE_FN

      integer    MPI_COMM_NULL_COPY_FN
      integer    MPI_COMM_DUP_FN
      integer    MPI_COMM_NULL_DELETE_FN
      common /MPIPRIVCOMATRCOPY/ MPI_COMM_NULL_COPY_FN, MPI_COMM_DUP_FN
      common /MPIPRIVCOMATRDEL/  MPI_COMM_NULL_DELETE_FN

! --- Empty group
      integer    MPI_GROUP_EMPTY
      parameter (MPI_GROUP_EMPTY = 0)

! --- Topologies
      integer    MPI_GRAPH
      parameter (MPI_GRAPH = 100)
      integer    MPI_CART
      parameter (MPI_CART = 101)

! --- MPI-2 Assertion Constants
      integer    MPI_MODE_NOCHECK
      parameter (MPI_MODE_NOCHECK      = 1)
      integer    MPI_MODE_NOSTORE
      parameter (MPI_MODE_NOSTORE      = 2)
      integer    MPI_MODE_NOPUT
      parameter (MPI_MODE_NOPUT        = 4)
      integer    MPI_MODE_NOPRECEDE
      parameter (MPI_MODE_NOPRECEDE    = 8)
      integer    MPI_MODE_NOSUCCEED
      parameter (MPI_MODE_NOSUCCEED    = 16)

! --- MPI-2 Locktype Constants
      integer    MPI_LOCK_EXCLUSIVE
      parameter (MPI_LOCK_EXCLUSIVE    = 1)
      integer    MPI_LOCK_SHARED
      parameter (MPI_LOCK_SHARED       = 2)

! --- Timing Functions
! --- (should be DOUBLE PRECISION -- not on Cray T3D/T3E)
      external MPI_WTICK
      real*8   MPI_WTICK
      external MPI_WTIME
      real*8   MPI_WTIME
      external PMPI_WTICK
      real*8   PMPI_WTICK
      external PMPI_WTIME
      real*8   PMPI_WTIME

      integer MPI_OFFSET_KIND
      parameter (MPI_OFFSET_KIND        = 8)

      integer MPI_ADDRESS_KIND
      parameter (MPI_ADDRESS_KIND       = 8)


! info parameters

      INTEGER(KIND=MPI_ADDRESS_KIND), PARAMETER:: MPI_INFO_NULL=0
      INTEGER, PARAMETER :: MPI_MAX_INFO_KEY=255
      INTEGER, PARAMETER :: MPI_MAX_INFO_VAL=1024


! MPI-2 I/O definitions

      include 'mpiof.h'
