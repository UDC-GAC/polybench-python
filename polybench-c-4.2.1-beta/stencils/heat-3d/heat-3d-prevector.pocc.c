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
/* heat-3d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "heat-3d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		 DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int i, j, k;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
        A[i][j][k] = B[i][j][k] = (DATA_TYPE) (i + j + (n-k))* 10 / (n);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n))

{
  int i, j, k;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++) {
         if ((i * n * n + j * n + k) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
         fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j][k]);
      }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_heat_3d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		      DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int t, i, j, k;

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
  register int c0, c2, c3, c1;
#pragma scop
if ((TSTEPS >= 1) && (_PB_N >= 3)) {
  for (c0 = 1; c0 <= TSTEPS; c0++) {
    {
      lb1 = ((2 * c0) + 1);
      ub1 = (((2 * c0) + _PB_N) + -2);
#pragma omp parallel for private(c3) firstprivate(lb1, ub1)
      for (c2 = lb1; c2 <= ub1; c2++) {
        {
#pragma ivdep
#pragma vector always
#pragma simd
          for (c3 = ((2 * c0) + 1); c3 <= (((2 * c0) + _PB_N) + -2); c3++) {
            {
              B[1][(-2 * c0) + c2][(-2 * c0) + c3] = SCALAR_VAL(0.125) * (A[1 +1][(-2 * c0) + c2][(-2 * c0) + c3] - SCALAR_VAL(2.0) * A[1][(-2 * c0) + c2][(-2 * c0) + c3] + A[1 -1][(-2 * c0) + c2][(-2 * c0) + c3]) + SCALAR_VAL(0.125) * (A[1][(-2 * c0) + c2+1][(-2 * c0) + c3] - SCALAR_VAL(2.0) * A[1][(-2 * c0) + c2][(-2 * c0) + c3] + A[1][(-2 * c0) + c2-1][(-2 * c0) + c3]) + SCALAR_VAL(0.125) * (A[1][(-2 * c0) + c2][(-2 * c0) + c3+1] - SCALAR_VAL(2.0) * A[1][(-2 * c0) + c2][(-2 * c0) + c3] + A[1][(-2 * c0) + c2][(-2 * c0) + c3-1]) + A[1][(-2 * c0) + c2][(-2 * c0) + c3];
            }
          }
        }
      }
      for (c1 = ((2 * c0) + 2); c1 <= (((2 * c0) + _PB_N) + -2); c1++) {
        {
#pragma ivdep
#pragma vector always
#pragma simd
          for (c3 = ((2 * c0) + 1); c3 <= (((2 * c0) + _PB_N) + -2); c3++) {
            {
              B[(-2 * c0) + c1][1][(-2 * c0) + c3] = SCALAR_VAL(0.125) * (A[(-2 * c0) + c1+1][1][(-2 * c0) + c3] - SCALAR_VAL(2.0) * A[(-2 * c0) + c1][1][(-2 * c0) + c3] + A[(-2 * c0) + c1-1][1][(-2 * c0) + c3]) + SCALAR_VAL(0.125) * (A[(-2 * c0) + c1][1 +1][(-2 * c0) + c3] - SCALAR_VAL(2.0) * A[(-2 * c0) + c1][1][(-2 * c0) + c3] + A[(-2 * c0) + c1][1 -1][(-2 * c0) + c3]) + SCALAR_VAL(0.125) * (A[(-2 * c0) + c1][1][(-2 * c0) + c3+1] - SCALAR_VAL(2.0) * A[(-2 * c0) + c1][1][(-2 * c0) + c3] + A[(-2 * c0) + c1][1][(-2 * c0) + c3-1]) + A[(-2 * c0) + c1][1][(-2 * c0) + c3];
            }
          }
          for (c2 = ((2 * c0) + 2); c2 <= (((2 * c0) + _PB_N) + -2); c2++) {
            {
              B[(-2 * c0) + c1][(-2 * c0) + c2][1] = SCALAR_VAL(0.125) * (A[(-2 * c0) + c1+1][(-2 * c0) + c2][1] - SCALAR_VAL(2.0) * A[(-2 * c0) + c1][(-2 * c0) + c2][1] + A[(-2 * c0) + c1-1][(-2 * c0) + c2][1]) + SCALAR_VAL(0.125) * (A[(-2 * c0) + c1][(-2 * c0) + c2+1][1] - SCALAR_VAL(2.0) * A[(-2 * c0) + c1][(-2 * c0) + c2][1] + A[(-2 * c0) + c1][(-2 * c0) + c2-1][1]) + SCALAR_VAL(0.125) * (A[(-2 * c0) + c1][(-2 * c0) + c2][1 +1] - SCALAR_VAL(2.0) * A[(-2 * c0) + c1][(-2 * c0) + c2][1] + A[(-2 * c0) + c1][(-2 * c0) + c2][1 -1]) + A[(-2 * c0) + c1][(-2 * c0) + c2][1];
#pragma ivdep
#pragma vector always
#pragma simd
              for (c3 = ((2 * c0) + 2); c3 <= (((2 * c0) + _PB_N) + -2); c3++) {
                {
                  A[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1] = SCALAR_VAL(0.125) * (B[((-2 * c0) + c1) + -1 +1][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1] - SCALAR_VAL(2.0) * B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1] + B[((-2 * c0) + c1) + -1 -1][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1]) + SCALAR_VAL(0.125) * (B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1 +1][((-2 * c0) + c3) + -1] - SCALAR_VAL(2.0) * B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1] + B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1 -1][((-2 * c0) + c3) + -1]) + SCALAR_VAL(0.125) * (B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1 +1] - SCALAR_VAL(2.0) * B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1] + B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1 -1]) + B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1];
                  B[(-2 * c0) + c1][(-2 * c0) + c2][(-2 * c0) + c3] = SCALAR_VAL(0.125) * (A[(-2 * c0) + c1+1][(-2 * c0) + c2][(-2 * c0) + c3] - SCALAR_VAL(2.0) * A[(-2 * c0) + c1][(-2 * c0) + c2][(-2 * c0) + c3] + A[(-2 * c0) + c1-1][(-2 * c0) + c2][(-2 * c0) + c3]) + SCALAR_VAL(0.125) * (A[(-2 * c0) + c1][(-2 * c0) + c2+1][(-2 * c0) + c3] - SCALAR_VAL(2.0) * A[(-2 * c0) + c1][(-2 * c0) + c2][(-2 * c0) + c3] + A[(-2 * c0) + c1][(-2 * c0) + c2-1][(-2 * c0) + c3]) + SCALAR_VAL(0.125) * (A[(-2 * c0) + c1][(-2 * c0) + c2][(-2 * c0) + c3+1] - SCALAR_VAL(2.0) * A[(-2 * c0) + c1][(-2 * c0) + c2][(-2 * c0) + c3] + A[(-2 * c0) + c1][(-2 * c0) + c2][(-2 * c0) + c3-1]) + A[(-2 * c0) + c1][(-2 * c0) + c2][(-2 * c0) + c3];
                }
              }
              A[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][_PB_N + -2] = SCALAR_VAL(0.125) * (B[((-2 * c0) + c1) + -1 +1][((-2 * c0) + c2) + -1][_PB_N + -2] - SCALAR_VAL(2.0) * B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][_PB_N + -2] + B[((-2 * c0) + c1) + -1 -1][((-2 * c0) + c2) + -1][_PB_N + -2]) + SCALAR_VAL(0.125) * (B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1 +1][_PB_N + -2] - SCALAR_VAL(2.0) * B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][_PB_N + -2] + B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1 -1][_PB_N + -2]) + SCALAR_VAL(0.125) * (B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][_PB_N + -2 +1] - SCALAR_VAL(2.0) * B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][_PB_N + -2] + B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][_PB_N + -2 -1]) + B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1][_PB_N + -2];
            }
          }
#pragma ivdep
#pragma vector always
#pragma simd
          for (c3 = ((2 * c0) + 2); c3 <= (((2 * c0) + _PB_N) + -1); c3++) {
            {
              A[((-2 * c0) + c1) + -1][_PB_N + -2][((-2 * c0) + c3) + -1] = SCALAR_VAL(0.125) * (B[((-2 * c0) + c1) + -1 +1][_PB_N + -2][((-2 * c0) + c3) + -1] - SCALAR_VAL(2.0) * B[((-2 * c0) + c1) + -1][_PB_N + -2][((-2 * c0) + c3) + -1] + B[((-2 * c0) + c1) + -1 -1][_PB_N + -2][((-2 * c0) + c3) + -1]) + SCALAR_VAL(0.125) * (B[((-2 * c0) + c1) + -1][_PB_N + -2 +1][((-2 * c0) + c3) + -1] - SCALAR_VAL(2.0) * B[((-2 * c0) + c1) + -1][_PB_N + -2][((-2 * c0) + c3) + -1] + B[((-2 * c0) + c1) + -1][_PB_N + -2 -1][((-2 * c0) + c3) + -1]) + SCALAR_VAL(0.125) * (B[((-2 * c0) + c1) + -1][_PB_N + -2][((-2 * c0) + c3) + -1 +1] - SCALAR_VAL(2.0) * B[((-2 * c0) + c1) + -1][_PB_N + -2][((-2 * c0) + c3) + -1] + B[((-2 * c0) + c1) + -1][_PB_N + -2][((-2 * c0) + c3) + -1 -1]) + B[((-2 * c0) + c1) + -1][_PB_N + -2][((-2 * c0) + c3) + -1];
            }
          }
        }
      }
      lb1 = ((2 * c0) + 2);
      ub1 = (((2 * c0) + _PB_N) + -1);
#pragma omp parallel for private(c3) firstprivate(lb1, ub1)
      for (c2 = lb1; c2 <= ub1; c2++) {
        {
#pragma ivdep
#pragma vector always
#pragma simd
          for (c3 = ((2 * c0) + 2); c3 <= (((2 * c0) + _PB_N) + -1); c3++) {
            {
              A[_PB_N + -2][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1] = SCALAR_VAL(0.125) * (B[_PB_N + -2 +1][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1] - SCALAR_VAL(2.0) * B[_PB_N + -2][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1] + B[_PB_N + -2 -1][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1]) + SCALAR_VAL(0.125) * (B[_PB_N + -2][((-2 * c0) + c2) + -1 +1][((-2 * c0) + c3) + -1] - SCALAR_VAL(2.0) * B[_PB_N + -2][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1] + B[_PB_N + -2][((-2 * c0) + c2) + -1 -1][((-2 * c0) + c3) + -1]) + SCALAR_VAL(0.125) * (B[_PB_N + -2][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1 +1] - SCALAR_VAL(2.0) * B[_PB_N + -2][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1] + B[_PB_N + -2][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1 -1]) + B[_PB_N + -2][((-2 * c0) + c2) + -1][((-2 * c0) + c3) + -1];
            }
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
  int n = N;
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_3D_ARRAY_DECL(A, DATA_TYPE, N, N, N, n, n, n);
  POLYBENCH_3D_ARRAY_DECL(B, DATA_TYPE, N, N, N, n, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_heat_3d (tsteps, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}
