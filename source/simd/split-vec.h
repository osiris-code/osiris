/*****************************************************************************************

The routines vsplit?D split particle trajectories into segments lying in the same cell for
use in the current deposition routines.

The current splitters expect that the particle paths stored in the vpbuf array are relative to the
cell where the motion started. This means that for particles crossing a cell boundary the values
will be outside the correct range. This is done to improve numerical accuracy since the total 
distance travelled is used as the denominator in a division. For small distances occuring next to
cell crossings, numerical roundoff could lead to divisions by 0 otherwise.

*****************************************************************************************/

#ifndef _SPLIT_VEC_H
#define _SPLIT_VEC_H

// If REAL has not been defined, do it automatically from the 
//PRECISION_SINGLE | PRECISION_DOUBLE macros

#ifndef REAL

#if defined( PRECISION_SINGLE )
#define REAL float
#elif defined( PRECISION_DOUBLE )
#define REAL double
#else
#error if REAL is not defined, PRECISION_SINGLE or PRECISION_DOUBLE must be defined.
#endif 

#endif

#ifndef VEC_WIDTH
#error VEC_WIDTH must be defined
#endif

/*****************************************************************************************
vsplit2D

Splits VEC_WIDTH particle trajectories for current deposition and stores segment on a special virtual
particle buffer which is optimized for the vector current deposition routines

This version uses a buffer which is a structure of arrays (SOA)

*****************************************************************************************/

void vsplit1D( t_split_buf1D* const vp, const int dix[] )
{

	unsigned int k, np;

	np = VEC_WIDTH;

	REAL* __restrict__ x0 = &(vp -> x0[ vp->np ]);
	REAL* __restrict__ x1 = &(vp -> x1[ vp->np ]);
	REAL* __restrict__  q = &(vp ->  q[ vp->np ]);
	REAL* __restrict__ vy = &(vp -> vy[ vp->np ]);
	REAL* __restrict__ vz = &(vp -> vz[ vp->np ]);

	int* __restrict__ ix = &(vp -> ix[ vp->np ]);

	for( k = 0 ; k < VEC_WIDTH; k++ ) {
		REAL delta, xint;

	    // Action is only required if dix != 0

		if ( dix[k] ) {

#if defined( PRECISION_SINGLE )
			xint  = 0.5f * dix[k];
#else
			xint  = 0.5 * dix[k];
#endif

			delta = ( xint - x0[k] ) / (x1[k] - x0[k]);

			// add extra particle (2nd segment)
			x0[np] = -xint; 	x1[np] = x1[k] - dix[k];
			ix[np] = ix[k] + dix[k]; 
			q[np]  = q[k];  
			vy[np] = vy[k] * (1-delta);
			vz[np] = vz[k] * (1-delta);
			np++;

			// correct existing particle (1st segment)
			x1[k] = xint;
			vy[k] *= delta;
			vz[k] *= delta;
		}  

	}

	vp -> np += np;
}

/*****************************************************************************************
vsplit2D

Splits VEC_WIDTH particle trajectories for current deposition and stores segment on a special virtual
particle buffer which is optimized for the vector current deposition routines

This version uses a buffer which is a structure of arrays (SOA)

*****************************************************************************************/

void vsplit2D( t_split_buf2D* const vp, const unsigned int cross[], const int dix[], const int diy[] )
{
 
  unsigned int k, np;

  np = VEC_WIDTH;
  
  REAL* __restrict__ x0 = &(vp -> x0[ vp->np ]);
  REAL* __restrict__ x1 = &(vp -> x1[ vp->np ]);
  REAL* __restrict__ y0 = &(vp -> y0[ vp->np ]);
  REAL* __restrict__ y1 = &(vp -> y1[ vp->np ]);
  REAL* __restrict__  q = &(vp ->  q[ vp->np ]);
  REAL* __restrict__ vz = &(vp -> vz[ vp->np ]);
  
  int* __restrict__ ix = &(vp -> ix[ vp->np ]);
  int* __restrict__ iy = &(vp -> iy[ vp->np ]);
  
  for( k = 0 ; k < VEC_WIDTH; k++ ) {
     REAL delta, xint, yint, vzint, xint2, yint2, vzint2;
          
     switch ( cross[k] ) {
       case (0) :  // no split

         // no action required
		 break;
         
       case (1) :  //  x cross only

		 xint  = 0.5 * dix[k];
		 delta = ( xint - x0[k] ) / (x1[k] - x0[k]);
		 yint  = y0[k] + (y1[k] - y0[k]) * delta;

         // add extra particle (2nd segment)
         x0[np] = -xint; 	x1[np] = x1[k] - dix[k];
         y0[np] = yint; 		y1[np] = y1[k];
         ix[np] = ix[k] + dix[k]; 
         iy[np] = iy[k];
         q[np]  = q[k];  
         vz[np] = vz[k] * (1-delta);
         np++;
         
         // correct existing particle (1st segment)
         x1[k] = xint;
         y1[k] = yint;
         vz[k] *= delta;

		 break;
          
       case (2) :  // y cross only

		 yint  = 0.5 * ( diy[k] );
		 delta = ( yint - y0[k] ) / (y1[k]-y0[k]);
		 xint  = x0[k] + (x1[k] - x0[k]) * delta;

         // add extra particle (2nd segment)
         x0[np] =  xint; 	x1[np] = x1[k];
         y0[np] = -yint; 	y1[np] = y1[k]  - diy[k];
         ix[np] = ix[k]; 
         iy[np] = iy[k] + diy[k];
         q[np]  = q[k];  
         vz[np] = vz[k] * (1-delta);
         np++;

         // correct existing particle (1st segment)         	 
		 x1[k] = xint;
		 y1[k] = yint;
		 vz[k] *= delta;
		 
		 break;
       
       case (3) : // x-y cross
		 // split in x direction first
		 xint  = 0.5 * dix[k];
		 delta = ( xint - x0[k] ) / (x1[k] - x0[k]);
		 yint  = y0[k] + (y1[k] - y0[k]) * delta;

         if ( ( yint >= -0.5 ) && (yint < 0.5) ) {

           // no y cross on 1st vp
           vzint = vz[k] * delta;

           // y split 2nd vp
           vzint2 = vz[k] * (1-delta);
           yint2 = 0.5 * diy[k];
		   delta = ( yint2 - yint ) / (y1[k] - yint);
		   xint2 = -xint + ( x1[k] - xint ) * delta;

           // add extra particle (1st segment y-split)
		   x0[np] =  -xint; 	x1[np] = xint2;
		   y0[np] = yint; 	y1[np] = yint2;
		   ix[np] = ix[k] + dix[k]; 
		   iy[np] = iy[k];
		   q[np]  = q[k];  
		   vz[np] = vzint2 * delta;
		   np++;


           // add extra particle (2nd segment y-split)
		   x0[np] =  xint2; 	x1[np] = x1[k] - dix[k];
		   y0[np] = -yint2; 	y1[np] = y1[k] - diy[k];
		   ix[np] = ix[k] + dix[k]; 
		   iy[np] = iy[k] + diy[k];
		   q[np]  = q[k];  
		   vz[np] = vzint2 * (1-delta);
		   np++;

		   // correct existing particle (1st segment)  
		   x1[k] = xint;
		   y1[k] = yint;
		   vz[k] = vzint;

         } else {
         
           vzint = vz[k]*delta;
		   // y split 1st vp
		   yint2 = 0.5 * diy[k];
		   delta = ( yint2 - y0[k] ) / ( yint - y0[k]);
		   xint2 = x0[k] + (xint - x0[k]) * delta;

           // add extra particle (2nd segment y-split)
		   x0[np] =  xint2; x1[np] = xint;
		   y0[np] = -yint2; y1[np] = yint - diy[k];
		   ix[np] = ix[k]; 
		   iy[np] = iy[k] + diy[k];
		   q[np]  = q[k];  
		   vz[np] = vzint * (1-delta);
		   np++;
 
		   // no y cross on last particle
		   x0[np] = -xint;          x1[np] = x1[k] - dix[k];
		   y0[np] = yint - diy[k];  y1[np] = y1[k] - diy[k];
		   ix[np] = ix[k] + dix[k];
		   iy[np] = iy[k] + diy[k];
		   q[np]  = q[k];
		   vz[np] = vz[k] - vzint;
		   np++;

		   // correct existing particle (1st segment y-split)  
		   x1[k] = xint2;
		   y1[k] = yint2;
		   vz[k] = vzint * delta;

         }
		 break;
     }
     
  }

  vp -> np += np;
}

/*****************************************************************************************
vsplit3D

Splits VEC_WIDTH particle trajectories for current deposition and stores segment on a special virtual
particle buffer which is optimized for the vector current deposition routines

This version uses a buffer which is a structure of arrays (SOA)

*****************************************************************************************/

void vsplit3D( t_split_buf3D* const vp, const unsigned int cross[], 
               const int dix[], const int diy[], const int diz[] )
{

  unsigned int k, np;

  np = VEC_WIDTH;
  
  REAL* __restrict__ x0 = &(vp -> x0[ vp->np ]);
  REAL* __restrict__ x1 = &(vp -> x1[ vp->np ]);
  REAL* __restrict__ y0 = &(vp -> y0[ vp->np ]);
  REAL* __restrict__ y1 = &(vp -> y1[ vp->np ]);
  REAL* __restrict__ z0 = &(vp -> z0[ vp->np ]);
  REAL* __restrict__ z1 = &(vp -> z1[ vp->np ]);
  REAL* __restrict__  q = &(vp ->  q[ vp->np ]);
  
  int* __restrict__ ix = &(vp -> ix[ vp->np ]);
  int* __restrict__ iy = &(vp -> iy[ vp->np ]);
  int* __restrict__ iz = &(vp -> iz[ vp->np ]);

  for( k = 0; k < VEC_WIDTH; k++ ) {
     REAL xint, yint, zint;
     REAL xint2, yint2, zint2;
     REAL delta;
                     
     switch ( cross[k] ) {
       case (0) :  // no split

         // no action required
		 break;
         
       case (1) :  //  x cross only
                  
		 xint  = 0.5 * dix[k];
		 delta = ( xint - x0[k] ) / (x1[k] - x0[k]);
		 yint  = y0[k] + (y1[k] - y0[k]) * delta;
		 zint  = z0[k] + (z1[k] - z0[k]) * delta;

         // add extra particle (2nd segment)
         x0[np] = -xint; 		x1[np] = x1[k] - dix[k];
         y0[np] =  yint; 		y1[np] = y1[k];
         z0[np] =  zint; 		z1[np] = z1[k];
         ix[np] = ix[k] + dix[k];
         iy[np] = iy[k];
         iz[np] = iz[k];
         q[np]  = q[k];
         np ++;

         // correct existing particle (1st segment)
         x1[k] = xint; 
         y1[k] = yint; 
         z1[k] = zint;

		 break;
          
       case (2) :  // y cross only
                  
		 yint  = 0.5 * diy[k];
		 delta = ( yint - y0[k] ) / (y1[k] - y0[k]);
		 xint  = x0[k] + (x1[k] - x0[k]) * delta;
		 zint  = z0[k] + (z1[k] - z0[k]) * delta;

         // add extra particle (2nd segment)
         x0[np] =  xint; 		x1[np] = x1[k];
         y0[np] = -yint; 		y1[np] = y1[k] - diy[k];
         z0[np] =  zint; 		z1[np] = z1[k];
         ix[np] = ix[k];
         iy[np] = iy[k] + diy[k];
         iz[np] = iz[k];
         q[np]  = q[k];
         np ++;

         // correct existing particle (1st segment)
         x1[k] = xint; 
         y1[k] = yint; 
         z1[k] = zint;

		 break;
       
       case (3) : // x-y cross
       
		 // split in x direction first
		 xint  = 0.5 * dix[k];
		 delta = ( xint - x0[k] ) / (x1[k] - x0[k]);

		 // note that this is different from case(1) because of change of cell in y direction
		 yint  = y0[k] + (y1[k] - y0[k]) * delta;
		 zint  = z0[k] + (z1[k] - z0[k]) * delta;

         if ( ( yint >= -0.5 ) && (yint < 0.5) ) {
                      
           // y split 2nd vp
           yint2 = 0.5 * diy[k];
		   delta = ( yint2 - yint ) / (y1[k] - yint);
		   xint2 = -xint + ( x1[k] - xint ) * delta;
		   zint2 =  zint + ( z1[k] - zint ) * delta;
		   
		   x0[np] = -xint; 		x1[np] = xint2;
		   y0[np] =  yint; 		y1[np] = yint2;
		   z0[np] =  zint; 		z1[np] = zint2;
		   ix[np] = ix[k] + dix[k];
		   iy[np] = iy[k];
		   iz[np] = iz[k];
		   q[np]  = q[k];
		   np ++;

		   x0[np] =  xint2; 		x1[np] = x1[k] - dix[k];
		   y0[np] = -yint2; 		y1[np] = y1[k] - diy[k];
		   z0[np] =  zint2; 		z1[np] = z1[k];
		   ix[np] = ix[k] + dix[k];
		   iy[np] = iy[k] + diy[k];
		   iz[np] = iz[k];
		   q[np]  = q[k];
		   np ++;

           // no y cross on 1st vp
           x1[k] = xint; y1[k] = yint; z1[k] = zint;
           
         } else {
         
		   // y split 1st vp
		   yint2 = 0.5 * (diy[k]);
		   delta = ( yint2 - y0[k] ) / ( yint - y0[k]);
		   xint2 = x0[k] + (xint - x0[k]) * delta;
		   zint2 = z0[k] + (zint - z0[k]) * delta;

		   x0[np] =  xint2; 		x1[np] = xint;
		   y0[np] = -yint2; 		y1[np] = yint - diy[k];
		   z0[np] =  zint2; 		z1[np] = zint;
		   ix[np] = ix[k];
		   iy[np] = iy[k] + diy[k];
		   iz[np] = iz[k];
		   q[np]  = q[k];
		   np ++;
					  
		   // no y cross on 2nd vp
		   x0[np] =       -xint;	x1[np] = x1[k] - dix[k];
		   y0[np] = yint-diy[k]; y1[np] = y1[k] - diy[k];
		   z0[np] =        zint;	z1[np] = z1[k];
		   ix[np] = ix[k] + dix[k];
		   iy[np] = iy[k] + diy[k];
		   iz[np] = iz[k];
		   q[np]  = q[k];
		   np ++;

           // Correct 1st vp
           x1[k] = xint2; y1[k] = yint2; z1[k] = zint2;

         }
  
		 break;

       case (4) :  // z cross only

		 zint  = 0.5 * diz[k];
		 delta = ( zint - z0[k] ) / (z1[k]-z0[k]);
		 xint  = x0[k] + (x1[k] - x0[k]) * delta;
		 yint  = y0[k] + (y1[k] - y0[k]) * delta;

         // add extra particle (2nd segment)
         x0[np] =  xint; 		x1[np] = x1[k];
         y0[np] =  yint; 		y1[np] = y1[k];
         z0[np] = -zint; 		z1[np] = z1[k] - diz[k];
         ix[np] = ix[k];
         iy[np] = iy[k];
         iz[np] = iz[k] + diz[k];
         q[np]  = q[k];
         np ++;

         // correct existing particle (1st segment)
         x1[k] = xint; y1[k] = yint; z1[k] = zint;
	 
		 break;

       case (5) :  // x-z cross
       
		 // split in x direction first
		 xint  = 0.5*( dix[k] );
		 delta = ( xint - x0[k] ) / (x1[k] - x0[k]);

		 // note that this is different from case(1) because of change of cell in z direction
		 yint  = y0[k] + (y1[k] - y0[k]) * delta; // no y cross
		 zint  = z0[k] + (z1[k] - z0[k]) * delta; 

         if ( ( zint >= -0.5 ) && (zint < 0.5) ) {
                      
           // z split 2nd vp
           zint2 = 0.5 * (diz[k]);
		   delta = ( zint2 - zint ) / (z1[k] - zint);
		   xint2 = -xint + ( x1[k] - xint ) * delta;
		   yint2 =  yint + ( y1[k] - yint ) * delta;
	
		   x0[np] = -xint; 		x1[np] = xint2;
		   y0[np] =  yint; 		y1[np] = yint2;
		   z0[np] =  zint; 		z1[np] = zint2;
		   ix[np] = ix[k] + dix[k];
		   iy[np] = iy[k];
		   iz[np] = iz[k];
		   q[np]  = q[k];
		   np ++;

		   x0[np] =  xint2; 		x1[np] = x1[k] - dix[k];
		   y0[np] =  yint2; 		y1[np] = y1[k];
		   z0[np] = -zint2; 		z1[np] = z1[k] - diz[k];
		   ix[np] = ix[k] + dix[k];
		   iy[np] = iy[k];
		   iz[np] = iz[k] + diz[k];
		   q[np]  = q[k];
		   np ++;

           // no z cross on 1st vp
		   x1[k] = xint;
		   y1[k] = yint;
		   z1[k] = zint;

         } else {
         
		   // z split 1st vp
		   zint2 = 0.5 * (diz[k]);
		   delta = ( zint2 - z0[k] ) / ( zint - z0[k]);
		   xint2 = x0[k] + (xint - x0[k]) * delta;
		   yint2 = y0[k] + (yint - y0[k]) * delta;

		   x0[np] =  xint2; 		x1[np] = xint;
		   y0[np] =  yint2; 		y1[np] = yint;
		   z0[np] = -zint2; 		z1[np] = zint - diz[k];
		   ix[np] = ix[k];
		   iy[np] = iy[k];
		   iz[np] = iz[k] + diz[k];
		   q[np]  = q[k];
		   np ++;
           					  
		   // no y cross on 2nd vp
		   x0[np] = -xint; 		x1[np] = x1[k] - dix[k];
		   y0[np] =  yint; 		y1[np] = y1[k];
		   z0[np] = zint-diz[k];	z1[np] = z1[k] - diz[k];
		   ix[np] = ix[k] + dix[k];
		   iy[np] = iy[k];
		   iz[np] = iz[k] + diz[k];
		   q[np]  = q[k];
		   np ++;

           // correct existing particle
		   x1[k] = xint2;
		   y1[k] = yint2;
		   z1[k] = zint2;

         }
		 break;       

       case (6) :  // y-z cross
       
		 // split in x direction first
		 yint  = 0.5*( diy[k] );
		 delta = ( yint - y0[k] ) / (y1[k] - y0[k]);

		 xint  = x0[k] + (x1[k] - x0[k]) * delta; // no x cross
		 zint  = z0[k] + (z1[k] - z0[k]) * delta; 

         if ( ( zint >= -0.5 ) && (zint < 0.5) ) {
                      
           // z split 2nd vp
           zint2 = 0.5 * diz[k];
		   delta = ( zint2 - zint ) / (z1[k] - zint);
		   xint2 =   xint + ( x1[k] - xint ) * delta;
		   yint2 =  -yint + ( y1[k] - yint ) * delta;
		   
		   x0[np] =  xint; 		x1[np] = xint2;
		   y0[np] = -yint; 		y1[np] = yint2;
		   z0[np] =  zint; 		z1[np] = zint2;
		   ix[np] = ix[k];
		   iy[np] = iy[k] + diy[k];
		   iz[np] = iz[k];
		   q[np]  = q[k];
		   np ++;

		   x0[np] =  xint2; 		x1[np] = x1[k];
		   y0[np] =  yint2; 		y1[np] = y1[k] - diy[k];
		   z0[np] = -zint2;	 	z1[np] = z1[k] - diz[k];
		   ix[np] = ix[k];
		   iy[np] = iy[k] + diy[k];
		   iz[np] = iz[k] + diz[k];
		   q[np]  = q[k];
		   np ++;

           // no z cross on 1st vp
		   x1[k] = xint; 
		   y1[k] = yint;
		   z1[k] = zint;
           
         } else {
         
		   // z split 1st vp
		   zint2 = 0.5 * diz[k];
		   delta = ( zint2 - z0[k] ) / ( zint - z0[k]);
		   xint2 = x0[k] + (xint - x0[k]) * delta;
		   yint2 = y0[k] + (yint - y0[k]) * delta;

		   x0[np] =  xint2; 		x1[np] = xint;
		   y0[np] =  yint2; 		y1[np] = yint;
		   z0[np] = -zint2; 		z1[np] = zint - diz[k];
		   ix[np] = ix[k];
		   iy[np] = iy[k];
		   iz[np] = iz[k] + diz[k];
		   q[np]  = q[k];
		   np ++;
					  
		   // no y cross on 2nd vp
		   x0[np] =  xint; 		x1[np] = x1[k];
		   y0[np] = -yint; 		y1[np] = y1[k] - diy[k];
		   z0[np] = zint-diz[k];	z1[np] = z1[k] - diz[k];
		   ix[np] = ix[k];
		   iy[np] = iy[k] + diy[k];
		   iz[np] = iz[k] + diz[k];
		   q[np]  = q[k];
		   np ++;

           // correct existing particle
		   x1[k] = xint2;
		   y1[k] = yint2;
		   z1[k] = zint2;

         }
		 break;       

       case (7) : // x-y-z cross
       {
        
          // do an x,y cross first 
		  xint  = 0.5 * dix[k];
		  delta = ( xint - x0[k] ) / (x1[k] - x0[k]);
		  yint  = y0[k] + (y1[k] - y0[k]) * delta;
		  zint  = z0[k] + (z1[k] - z0[k]) * delta; 

         if ( ( yint >= -0.5 ) && (yint < 0.5) ) {
                      
           // y split 2nd vp
           yint2 = 0.5 * (diy[k]);
		   delta = ( yint2 - yint ) / (y1[k] - yint);
		   xint2 = -xint + ( x1[k] - xint ) * delta;
		   zint2 =  zint + ( z1[k] - zint ) * delta;
		   
		   x0[np] = -xint; 		x1[np] = xint2;
		   y0[np] =  yint; 		y1[np] = yint2;
		   z0[np] =  zint; 		z1[np] = zint2;
		   ix[np] = ix[k] + dix[k];
		   iy[np] = iy[k];
		   iz[np] = iz[k];
		   q[np]  = q[k];
		   np ++;

		   x0[np] =  xint2; 		x1[np] = x1[k] - dix[k];
		   y0[np] = -yint2; 		y1[np] = y1[k] - diy[k];
		   z0[np] =  zint2; 		z1[np] = z1[k];
		   ix[np] = ix[k] + dix[k];
		   iy[np] = iy[k] + diy[k];
		   iz[np] = iz[k];
		   q[np]  = q[k];
		   np ++;

           // no y cross on 1st vp
           x1[k] = xint; 
           y1[k] = yint;
           z1[k] = zint;
           
         } else {
         
		   // y split 1st vp
		   yint2 = 0.5 * (diy[k]);
		   delta = ( yint2 - y0[k] ) / ( yint - y0[k]);
		   xint2 = x0[k] + (xint - x0[k]) * delta;
		   zint2 = z0[k] + (zint - z0[k]) * delta;

		   x0[np] =  xint2; 		x1[np] = xint;
		   y0[np] = -yint2; 		y1[np] = yint - diy[k];
		   z0[np] =  zint2; 		z1[np] = zint;
		   ix[np] = ix[k];
		   iy[np] = iy[k] + diy[k];
		   iz[np] = iz[k];
		   q[np]  = q[k];
		   np ++;
					  
		   // no y cross on 2nd vp
		   x0[np] =       -xint;	x1[np] = x1[k] - dix[k];
		   y0[np] = yint-diy[k]; y1[np] = y1[k] - diy[k];
		   z0[np] =        zint;	z1[np] = z1[k];
		   ix[np] = ix[k] + dix[k];
		   iy[np] = iy[k] + diy[k];
		   iz[np] = iz[k];
		   q[np]  = q[k];
		   np ++;

           // Correct 1st vp
           x1[k] = xint2; 
           y1[k] = yint2; 
           z1[k] = zint2;

         }
         
 		  // one of the 3 vp requires an additional z split
          // int vpidx[3] = {k, np-2, np-1};
          int vpidx[3]; vpidx[0] = k; vpidx[1] = np-2; vpidx[2] = np-1;
          
          int i, j, idx;

		  zint = 0.5 * (diz[k]);
		  for ( i = 0; i < 3; i ++ ) {
			idx = vpidx[i];
			if ( (z1[idx] < -0.5) || (z1[idx] >= +0.5) ) {
			   
			   // z1 and z0 are already indexed to the same cell
			   delta = ( zint - z0[idx]) / ( z1[idx]-z0[idx] );
			   xint  = x0[idx] + (x1[idx] - x0[idx]) * delta;
			   yint  = y0[idx] + (y1[idx] - y0[idx]) * delta;
	  
			   // store new vp
			   x0[np] =  xint; 		x1[np] = x1[idx];
			   y0[np] =  yint; 		y1[np] = y1[idx];
			   z0[np] = -zint; 		z1[np] = z1[idx] - diz[k];
			   ix[np] = ix[idx];
			   iy[np] = iy[idx];
			   iz[np] = iz[idx] + diz[k];
			   q[np]  = q[idx];
			   np ++;
			   		   
		       // correct old vp
			   x1[idx] = xint;
			   y1[idx] = yint;
			   z1[idx] = zint; 
			   
			   for ( j = i+1; j < 3; j++ ) {
			     idx = vpidx[j];
			     z0[idx] -= diz[k];
			     z1[idx] -= diz[k];
			     iz[idx] += diz[k];
			   }
			   
			   // No further z splits
			   break;			   
			}
		  }
		          
       }

     }
     
  }

  vp -> np += np;

}

#endif
