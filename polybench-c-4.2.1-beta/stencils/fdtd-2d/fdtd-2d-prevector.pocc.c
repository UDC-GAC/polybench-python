#include <math.h>
/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* fdtd-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "fdtd-2d.h"


/* Array initialization. */
static
void init_array (int tmax,
		 int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int i, j;

  for (i = 0; i < tmax; i++)
    _fict_[i] = (DATA_TYPE) i;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      {
	ex[i][j] = ((DATA_TYPE) i*(j+1)) / nx;
	ey[i][j] = ((DATA_TYPE) i*(j+2)) / ny;
	hz[i][j] = ((DATA_TYPE) i*(j+3)) / nx;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("ex");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ex[i][j]);
    }
  POLYBENCH_DUMP_END("ex");
  POLYBENCH_DUMP_FINISH;

  POLYBENCH_DUMP_BEGIN("ey");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ey[i][j]);
    }
  POLYBENCH_DUMP_END("ey");

  POLYBENCH_DUMP_BEGIN("hz");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, hz[i][j]);
    }
  POLYBENCH_DUMP_END("hz");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_fdtd_2d(int tmax,
		    int nx,
		    int ny,
		    DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int t, i, j;

#ifdef ceild
# undef ceild
#endif
#ifdef floord
# undef floord
#endif
#ifdef max
# undef max
#endif
#ifdef min
# undef min
#endif
#define ceild(n,d) (((n) < 0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define floord(x,y) (((x) < 0)? -((-(x)+(y)-1)/(y)) : (x)/(y))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))
/* Copyright (C) 1991-2020 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses Unicode 10.0.0.  Version 10.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2017, fifth edition, plus
   the following additions from Amendment 1 to the fifth edition:
   - 56 emoji characters
   - 285 hentaigana
   - 3 additional Zanabazar Square characters */
  register int lbv, ubv, lb, ub, lb1, ub1, lb2, ub2;
  register int c0, c2, c1;
#pragma scop
if ((_PB_NY >= 1) && (_PB_TMAX >= 1)) {
  if ((_PB_NX >= 2) && (_PB_NY >= 2)) {
    for (c0 = 0; c0 <= (_PB_TMAX + -1); c0++) {
      {
        ey[0][0] = _fict_[c0];
#pragma ivdep
#pragma vector always
#pragma simd
        for (c2 = (c0 + 1); c2 <= ((c0 + _PB_NX) + -1); c2++) {
          {
            ey[(-1 * c0) + c2][0] = ey[(-1 * c0) + c2][0] - SCALAR_VAL(0.5)*(hz[(-1 * c0) + c2][0]-hz[(-1 * c0) + c2-1][0]);
          }
        }
        for (c1 = (c0 + 1); c1 <= ((c0 + _PB_NY) + -1); c1++) {
          {
            ex[0][(-1 * c0) + c1] = ex[0][(-1 * c0) + c1] - SCALAR_VAL(0.5)*(hz[0][(-1 * c0) + c1]-hz[0][(-1 * c0) + c1-1]);
            ey[0][(-1 * c0) + c1] = _fict_[c0];
            for (c2 = (c0 + 1); c2 <= ((c0 + _PB_NX) + -1); c2++) {
              {
                hz[((-1 * c0) + c2) + -1][((-1 * c0) + c1) + -1] = hz[((-1 * c0) + c2) + -1][((-1 * c0) + c1) + -1] - SCALAR_VAL(0.7)* (ex[((-1 * c0) + c2) + -1][((-1 * c0) + c1) + -1 +1] - ex[((-1 * c0) + c2) + -1][((-1 * c0) + c1) + -1] + ey[((-1 * c0) + c2) + -1 +1][((-1 * c0) + c1) + -1] - ey[((-1 * c0) + c2) + -1][((-1 * c0) + c1) + -1]);
                ey[(-1 * c0) + c2][(-1 * c0) + c1] = ey[(-1 * c0) + c2][(-1 * c0) + c1] - SCALAR_VAL(0.5)*(hz[(-1 * c0) + c2][(-1 * c0) + c1]-hz[(-1 * c0) + c2-1][(-1 * c0) + c1]);
                ex[(-1 * c0) + c2][(-1 * c0) + c1] = ex[(-1 * c0) + c2][(-1 * c0) + c1] - SCALAR_VAL(0.5)*(hz[(-1 * c0) + c2][(-1 * c0) + c1]-hz[(-1 * c0) + c2][(-1 * c0) + c1-1]);
              }
            }
          }
        }
      }
    }
  }
  if ((_PB_NX >= 2) && (_PB_NY == 1)) {
    for (c0 = 0; c0 <= (_PB_TMAX + -1); c0++) {
      {
        ey[0][0] = _fict_[c0];
#pragma ivdep
#pragma vector always
#pragma simd
        for (c2 = (c0 + 1); c2 <= ((c0 + _PB_NX) + -1); c2++) {
          {
            ey[(-1 * c0) + c2][0] = ey[(-1 * c0) + c2][0] - SCALAR_VAL(0.5)*(hz[(-1 * c0) + c2][0]-hz[(-1 * c0) + c2-1][0]);
          }
        }
      }
    }
  }
  if ((_PB_NX == 1) && (_PB_NY >= 2)) {
    for (c0 = 0; c0 <= (_PB_TMAX + -1); c0++) {
      {
        ey[0][0] = _fict_[c0];
#pragma ivdep
#pragma vector always
#pragma simd
        for (c1 = (c0 + 1); c1 <= ((c0 + _PB_NY) + -1); c1++) {
          {
            ex[0][(-1 * c0) + c1] = ex[0][(-1 * c0) + c1] - SCALAR_VAL(0.5)*(hz[0][(-1 * c0) + c1]-hz[0][(-1 * c0) + c1-1]);
            ey[0][(-1 * c0) + c1] = _fict_[c0];
          }
        }
      }
    }
  }
  if ((_PB_NX <= 1) && (_PB_NY == 1)) {
    for (c0 = 0; c0 <= (_PB_TMAX + -1); c0++) {
      {
        ey[0][0] = _fict_[c0];
      }
    }
  }
  if ((_PB_NX <= 0) && (_PB_NY >= 2)) {
    for (c0 = 0; c0 <= (_PB_TMAX + -1); c0++) {
      {
#pragma ivdep
#pragma vector always
#pragma simd
        for (c1 = c0; c1 <= ((c0 + _PB_NY) + -1); c1++) {
          {
            ey[0][(-1 * c0) + c1] = _fict_[c0];
          }
        }
      }
    }
  }
}
#pragma endscop
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int tmax = TMAX;
  int nx = NX;
  int ny = NY;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(ex,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(ey,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(hz,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_1D_ARRAY_DECL(_fict_,DATA_TYPE,TMAX,tmax);

  /* Initialize array(s). */
  init_array (tmax, nx, ny,
	      POLYBENCH_ARRAY(ex),
	      POLYBENCH_ARRAY(ey),
	      POLYBENCH_ARRAY(hz),
	      POLYBENCH_ARRAY(_fict_));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_fdtd_2d (tmax, nx, ny,
		  POLYBENCH_ARRAY(ex),
		  POLYBENCH_ARRAY(ey),
		  POLYBENCH_ARRAY(hz),
		  POLYBENCH_ARRAY(_fict_));


  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(nx, ny, POLYBENCH_ARRAY(ex),
				    POLYBENCH_ARRAY(ey),
				    POLYBENCH_ARRAY(hz)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(ex);
  POLYBENCH_FREE_ARRAY(ey);
  POLYBENCH_FREE_ARRAY(hz);
  POLYBENCH_FREE_ARRAY(_fict_);

  return 0;
}
