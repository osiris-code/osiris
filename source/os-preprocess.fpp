!
!  File: os-preprocess.fpp
!
!  A set of preprocessing macros for osiris, and the fpp/cpp preprocessor
!
!

#if defined (__IBMC__)||(__GNUC__)||(_WIN32)||(_WIN64)

! Macros for the IBM CPP/GNU preprocessor

#define FILE_LINE " [",__FILE__,":",__LINE__,"]"


! Assertion macro

#define ASSERT(exp) call asrt__(exp,#exp,__FILE__,__LINE__)

! Debug functions
#if defined (DEBUG_FILE) || defined (DEBUG_GLOBAL)

#define DEBUG(...) write(dbg_buf__,*) __VA_ARGS__; call dbg__(__FILE__,__LINE__)

#else

#define DEBUG(...)

#endif

#define DEBUG_EXP(exp) DEBUG(#exp," = ",exp)


! SCR hostname functions

#define SCR_ROOT(...) if (mpi_node()==0) print *,__VA_ARGS__
#define SCR_HOST(...) print *, "[",trim(hostname()),"] ",__VA_ARGS__
#define SCR_MPINODE(...) print *, "[",mpi_node(),"] ",__VA_ARGS__
#define SCR_EXP_MPINODE(exp) SCR_MPINODE(#exp," = ",exp)


! LOG functions

#define LOG(...) write(log_buf__,*) __VA_ARGS__;call log__
#define LOG_EXP(exp) LOG(#exp," = ",exp)

! ERROR functions

#define ERROR(...) write(err_buf__,*) __VA_ARGS__;call err__(__FILE__,__LINE__)
#define ERROR_EXP(exp) ERROR(#exp," = ",exp)

#define CHECK_ERROR( ierr, msg, code ) call check_error(ierr,msg,code,__FILE__,__LINE__)

! WARNING functions

#define WARNING(...) write(wrn_buf__,*) __VA_ARGS__;call wrn__(__FILE__,__LINE__)
#define WARNING_EXP(exp) WARNING(#exp," = ",exp)

! functions for restart io
#define restart_io_wr(var, handle, ierr) call restart_io_write(#var, var, handle, ierr)
#define restart_io_rd(var, handle, ierr) call restart_io_read(var, handle, ierr)

#else


#error Compiler not (yet) supported. Please use an older version of osiris
#error while your compiler is not supported.

#endif
