#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 4
#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT
//#define NOCOMPLEXFMTCHANGE


#define GGA_WITH_PLAN
#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)



// Calling convention:
//  c = comp_chirpcz(f,K,deltao,o)

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{

   mwSize L  = mxGetM(prhs[0]);
   mwSize W  = mxGetN(prhs[0]);
   mwSize K = (mwSize) mxGetScalar(prhs[1]);
   double deltao = mxGetScalar(prhs[2]);
   double o = mxGetScalar(prhs[3]);

   const LTFAT_TYPE* fPtr = (const LTFAT_TYPE*) mxGetData(prhs[0]);

   plhs[0] = ltfatCreateMatrix(K,W,LTFAT_MX_CLASSID,mxCOMPLEX);
   LTFAT_REAL _Complex* cPtr = (LTFAT_REAL _Complex*) mxGetPr(plhs[0]);

   LTFAT_NAME(chzt)(fPtr,L,W,K,deltao,o,cPtr);
   //LTFAT_NAME(chzt_fact)(fPtr,L,W,K,deltao,o,cPtr);

   return;
}
#endif
