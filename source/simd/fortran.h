/*

 This file includes the necessary definitions to ease the calling of C functions from
 Fortran.

*/

#ifndef _FORTRAN_H
#define _FORTRAN_H

/*

Define the function names according the the fortran compiler external
names. The default is to generate names with a single trailing underscore.

*/

#ifdef FORTRANNOUNDERSCORE
#define FNAME(A) A
#endif

#ifdef FORTRANDOUBLEUNDERSCORE
#ifdef FNAME
#error FNAME macro already defined, check your FORTRAN* defines.
#else
#define FNAME(A) A##__
#endif
#endif

#ifdef FORTRANSINGLEUNDERSCORE
#ifdef FNAME
#error FNAME macro already defined, check your FORTRAN* defines.
#else
#define FNAME(A) A##_
#endif
#endif

#ifndef FNAME
#error FNAME macro not defined, check your FORTRAN* defines.
#endif


#endif /* _FORTRAN_H */
