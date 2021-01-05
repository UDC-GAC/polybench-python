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
/* adi.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "adi.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(u,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	u[i][j] =  (DATA_TYPE)(i + n-j) / n;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(u,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("u");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, u[i][j]);
    }
  POLYBENCH_DUMP_END("u");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Based on a Fortran code fragment from Figure 5 of
 * "Automatic Data and Computation Decomposition on Distributed Memory Parallel Computers"
 * by Peizong Lee and Zvi Meir Kedem, TOPLAS, 2002
 */
static
void kernel_adi(int tsteps, int n,
		DATA_TYPE POLYBENCH_2D(u,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(v,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(p,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(q,N,N,n,n))
{
  int t, i, j;
  DATA_TYPE DX, DY, DT;
  DATA_TYPE B1, B2;
  DATA_TYPE mul1, mul2;
  DATA_TYPE a, b, c, d, e, f;

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
  register int c1, c17, c18;
#pragma scop
DY = SCALAR_VAL(1.0)/_PB_N;
DT = SCALAR_VAL(1.0)/_PB_TSTEPS;
B1 = SCALAR_VAL(2.0);
B2 = SCALAR_VAL(1.0);
mul2 = B2 * DT / (DY * DY);
e = SCALAR_VAL(1.0)+mul2;
d = -mul2 / SCALAR_VAL(2.0);
DX = SCALAR_VAL(1.0)/_PB_N;
f = d;
mul1 = B1 * DT / (DX * DX);
b = SCALAR_VAL(1.0)+mul1;
a = -mul1 / SCALAR_VAL(2.0);
c = a;
if (_PB_TSTEPS >= 3) {
  for (c1 = 1; c1 <= _PB_N; c1++) {
    {
#pragma ivdep
#pragma vector always
#pragma simd
      for (c17 = 1; c17 <= (_PB_TSTEPS + -2); c17++) {
        {
          v[_PB_N-1][c17] = SCALAR_VAL(1.0);
        }
      }
#pragma ivdep
#pragma vector always
#pragma simd
      for (c17 = 1; c17 <= (_PB_TSTEPS + -2); c17++) {
        {
          p[c17][0] = SCALAR_VAL(0.0);
        }
      }
      ub1 = (_PB_TSTEPS + -2);
#pragma omp parallel for private(c18) firstprivate(ub1)
      for (c17 = 1; c17 <= ub1; c17++) {
        {
          for (c18 = 1; c18 <= (_PB_TSTEPS + -2); c18++) {
            {
              p[c17][c18] = -c / (a*p[c17][c18-1]+b);
            }
          }
        }
      }
#pragma ivdep
#pragma vector always
#pragma simd
      for (c17 = 1; c17 <= (_PB_TSTEPS + -2); c17++) {
        {
          v[0][c17] = SCALAR_VAL(1.0);
        }
      }
#pragma ivdep
#pragma vector always
#pragma simd
      for (c17 = 1; c17 <= (_PB_TSTEPS + -2); c17++) {
        {
          q[c17][0] = v[0][c17];
        }
      }
      ub1 = (_PB_TSTEPS + -2);
#pragma omp parallel for private(c18) firstprivate(ub1)
      for (c17 = 1; c17 <= ub1; c17++) {
        {
          for (c18 = 1; c18 <= (_PB_TSTEPS + -2); c18++) {
            {
              q[c17][c18] = (-d*u[c18][c17-1]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*d)*u[c18][c17] - f*u[c18][c17+1]-a*q[c17][c18-1])/(a*p[c17][c18-1]+b);
            }
          }
        }
      }
#pragma ivdep
#pragma vector always
#pragma simd
      for (c17 = 1; c17 <= (_PB_TSTEPS + -2); c17++) {
        {
          u[c17][0] = SCALAR_VAL(1.0);
          p[c17][0] = SCALAR_VAL(0.0);
          u[c17][_PB_N-1] = SCALAR_VAL(1.0);
        }
      }
#pragma ivdep
#pragma vector always
#pragma simd
      for (c17 = 1; c17 <= (_PB_TSTEPS + -2); c17++) {
        {
          q[c17][0] = u[c17][0];
        }
      }
      ub1 = (_PB_TSTEPS + -2);
#pragma omp parallel for private(c18) firstprivate(ub1)
      for (c17 = 1; c17 <= ub1; c17++) {
        {
          for (c18 = 1; c18 <= (_PB_TSTEPS + -2); c18++) {
            {
              v[_PB_N-1-c18][c17] = p[c17][_PB_N-1-c18] * v[_PB_N-c18][c17] + q[c17][_PB_N-1-c18];
            }
          }
        }
      }
      ub1 = (_PB_TSTEPS + -2);
#pragma omp parallel for private(c18) firstprivate(ub1)
      for (c17 = 1; c17 <= ub1; c17++) {
        {
          for (c18 = 1; c18 <= (_PB_TSTEPS + -2); c18++) {
            {
              p[c17][c18] = -f / (d*p[c17][c18-1]+e);
            }
          }
        }
      }
      ub1 = (_PB_TSTEPS + -2);
#pragma omp parallel for private(c18) firstprivate(ub1)
      for (c17 = 1; c17 <= ub1; c17++) {
        {
          for (c18 = 1; c18 <= (_PB_TSTEPS + -2); c18++) {
            {
              q[c17][c18] = (-a*v[c17-1][c18]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*a)*v[c17][c18] - c*v[c17+1][c18]-d*q[c17][c18-1])/(d*p[c17][c18-1]+e);
            }
          }
        }
      }
      ub1 = (_PB_TSTEPS + -2);
#pragma omp parallel for private(c18) firstprivate(ub1)
      for (c17 = 1; c17 <= ub1; c17++) {
        {
          for (c18 = 1; c18 <= (_PB_TSTEPS + -2); c18++) {
            {
              u[c17][_PB_N-1-c18] = p[c17][_PB_N-1-c18] * u[c17][_PB_N-c18] + q[c17][_PB_N-1-c18];
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
  POLYBENCH_2D_ARRAY_DECL(u, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(v, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(p, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(q, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(u));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_adi (tsteps, n, POLYBENCH_ARRAY(u), POLYBENCH_ARRAY(v), POLYBENCH_ARRAY(p), POLYBENCH_ARRAY(q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(u)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(u);
  POLYBENCH_FREE_ARRAY(v);
  POLYBENCH_FREE_ARRAY(p);
  POLYBENCH_FREE_ARRAY(q);

  return 0;
}
