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
/* correlation.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "correlation.h"


/* Array initialization. */
static
void init_array (int m,
		 int n,
		 DATA_TYPE *float_n,
		 DATA_TYPE POLYBENCH_2D(data,N,M,n,m))
{
  int i, j;

  *float_n = (DATA_TYPE)N;

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      data[i][j] = (DATA_TYPE)(i*j)/M + i;

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m,
		 DATA_TYPE POLYBENCH_2D(corr,M,M,m,m))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("corr");
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, corr[i][j]);
    }
  POLYBENCH_DUMP_END("corr");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_correlation(int m, int n,
			DATA_TYPE float_n,
			DATA_TYPE POLYBENCH_2D(data,N,M,n,m),
			DATA_TYPE POLYBENCH_2D(corr,M,M,m,m),
			DATA_TYPE POLYBENCH_1D(mean,M,m),
			DATA_TYPE POLYBENCH_1D(stddev,M,m))
{
  int i, j, k;

  DATA_TYPE eps = SCALAR_VAL(0.1);


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
  register int c1, c2, c3;
#pragma scop
corr[_PB_M-1][_PB_M-1] = SCALAR_VAL(1.0);
ub1 = (_PB_M + -2);
#pragma omp parallel for private(c2) firstprivate(ub1)
for (c1 = 0; c1 <= ub1; c1++) {
  {
#pragma ivdep
#pragma vector always
#pragma simd
    for (c2 = (c1 + 1); c2 <= (_PB_M + -1); c2++) {
      {
        corr[c1][c2] = SCALAR_VAL(0.0);
      }
    }
  }
}
#pragma ivdep
#pragma vector always
#pragma simd
for (c1 = 0; c1 <= (_PB_M + -2); c1++) {
  {
    corr[c1][c1] = SCALAR_VAL(1.0);
    stddev[c1] = SCALAR_VAL(0.0);
    mean[c1] = SCALAR_VAL(0.0);
  }
}
if (_PB_M >= 1) {
  stddev[_PB_M + -1] = SCALAR_VAL(0.0);
  mean[_PB_M + -1] = SCALAR_VAL(0.0);
}
if (_PB_N >= 1) {
  {
    for (c2 = 0; c2 <= (_PB_N + -1); c2++) {
#pragma ivdep
#pragma vector always
#pragma simd
      for (c1 = 0; c1 <= (_PB_M + -1); c1++) {
        {
          mean[c1] += data[c2][c1];
        }
      }
    }
  }
}
#pragma ivdep
#pragma vector always
#pragma simd
for (c1 = 0; c1 <= (_PB_M + -1); c1++) {
  {
    mean[c1] /= float_n;
  }
}
if (_PB_N >= 1) {
  {
    for (c2 = 0; c2 <= (_PB_N + -1); c2++) {
#pragma ivdep
#pragma vector always
#pragma simd
      for (c1 = 0; c1 <= (_PB_M + -1); c1++) {
        {
          stddev[c1] += (data[c2][c1] - mean[c1]) * (data[c2][c1] - mean[c1]);
          data[c2][c1] -= mean[c1];
        }
      }
    }
  }
}
#pragma ivdep
#pragma vector always
#pragma simd
for (c1 = 0; c1 <= (_PB_M + -1); c1++) {
  {
    stddev[c1] /= float_n;
    stddev[c1] = SQRT_FUN(stddev[c1]);
    stddev[c1] = stddev[c1] <= eps ? SCALAR_VAL(1.0) : stddev[c1];
  }
}
if (_PB_M >= 1) {
  ub1 = (_PB_N + -1);
#pragma omp parallel for private(c2) firstprivate(ub1)
  for (c1 = 0; c1 <= ub1; c1++) {
    {
#pragma ivdep
#pragma vector always
#pragma simd
      for (c2 = 0; c2 <= (_PB_M + -1); c2++) {
        {
          data[c1][c2] /= SQRT_FUN(float_n) * stddev[c2];
        }
      }
    }
  }
}
if (_PB_N >= 1) {
  ub1 = (_PB_M + -2);
#pragma omp parallel for private(c3, c2) firstprivate(ub1)
  for (c1 = 0; c1 <= ub1; c1++) {
    {
      {
        for (c3 = 0; c3 <= (_PB_N + -1); c3++) {
#pragma ivdep
#pragma vector always
#pragma simd
          for (c2 = (c1 + 1); c2 <= (_PB_M + -1); c2++) {
            {
              corr[c1][c2] += (data[c3][c1] * data[c3][c2]);
            }
          }
        }
      }
    }
  }
}
ub1 = (_PB_M + -2);
#pragma omp parallel for private(c2) firstprivate(ub1)
for (c1 = 0; c1 <= ub1; c1++) {
  {
#pragma ivdep
#pragma vector always
#pragma simd
    for (c2 = (c1 + 1); c2 <= (_PB_M + -1); c2++) {
      {
        corr[c2][c1] = corr[c1][c2];
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
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE float_n;
  POLYBENCH_2D_ARRAY_DECL(data,DATA_TYPE,N,M,n,m);
  POLYBENCH_2D_ARRAY_DECL(corr,DATA_TYPE,M,M,m,m);
  POLYBENCH_1D_ARRAY_DECL(mean,DATA_TYPE,M,m);
  POLYBENCH_1D_ARRAY_DECL(stddev,DATA_TYPE,M,m);

  /* Initialize array(s). */
  init_array (m, n, &float_n, POLYBENCH_ARRAY(data));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_correlation (m, n, float_n,
		      POLYBENCH_ARRAY(data),
		      POLYBENCH_ARRAY(corr),
		      POLYBENCH_ARRAY(mean),
		      POLYBENCH_ARRAY(stddev));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, POLYBENCH_ARRAY(corr)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(data);
  POLYBENCH_FREE_ARRAY(corr);
  POLYBENCH_FREE_ARRAY(mean);
  POLYBENCH_FREE_ARRAY(stddev);

  return 0;
}
