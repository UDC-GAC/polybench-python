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
/* 2mm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "2mm.h"


/* Array initialization. */
static
void init_array(int ni, int nj, int nk, int nl,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) ((i*j+1) % ni) / ni;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+1) % nj) / nj;
  for (i = 0; i < nj; i++)
    for (j = 0; j < nl; j++)
      C[i][j] = (DATA_TYPE) ((i*(j+3)+1) % nl) / nl;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
      D[i][j] = (DATA_TYPE) (i*(j+2) % nk) / nk;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nl,
		 DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("D");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, D[i][j]);
    }
  POLYBENCH_DUMP_END("D");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_2mm(int ni, int nj, int nk, int nl,
		DATA_TYPE alpha,
		DATA_TYPE beta,
		DATA_TYPE POLYBENCH_2D(tmp,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j, k;

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
  register int c1, c2, c5;
#pragma scop
if (_PB_NI >= 1) {
  if ((_PB_NJ >= 1) && (_PB_NL >= 1)) {
    ub1 = (_PB_NI + -1);
#pragma omp parallel for private(c2) firstprivate(ub1)
    for (c1 = 0; c1 <= ub1; c1++) {
      {
#pragma ivdep
#pragma vector always
#pragma simd
        for (c2 = 0; c2 <= min((_PB_NJ + -1), (_PB_NL + -1)); c2++) {
          {
            D[c1][c2] *= beta;
            tmp[c1][c2] = SCALAR_VAL(0.0);
          }
        }
#pragma ivdep
#pragma vector always
#pragma simd
        for (c2 = _PB_NL; c2 <= (_PB_NJ + -1); c2++) {
          {
            tmp[c1][c2] = SCALAR_VAL(0.0);
          }
        }
#pragma ivdep
#pragma vector always
#pragma simd
        for (c2 = _PB_NJ; c2 <= (_PB_NL + -1); c2++) {
          {
            D[c1][c2] *= beta;
          }
        }
      }
    }
  }
  if ((_PB_NJ >= 1) && (_PB_NL <= 0)) {
    ub1 = (_PB_NI + -1);
#pragma omp parallel for private(c2) firstprivate(ub1)
    for (c1 = 0; c1 <= ub1; c1++) {
      {
#pragma ivdep
#pragma vector always
#pragma simd
        for (c2 = 0; c2 <= (_PB_NJ + -1); c2++) {
          {
            tmp[c1][c2] = SCALAR_VAL(0.0);
          }
        }
      }
    }
  }
  if ((_PB_NJ <= 0) && (_PB_NL >= 1)) {
    ub1 = (_PB_NI + -1);
#pragma omp parallel for private(c2) firstprivate(ub1)
    for (c1 = 0; c1 <= ub1; c1++) {
      {
#pragma ivdep
#pragma vector always
#pragma simd
        for (c2 = 0; c2 <= (_PB_NL + -1); c2++) {
          {
            D[c1][c2] *= beta;
          }
        }
      }
    }
  }
  if (((_PB_NJ >= 1) && (_PB_NK >= 1)) && (_PB_NL >= 1)) {
    ub1 = (_PB_NI + -1);
#pragma omp parallel for private(c2, c5) firstprivate(ub1)
    for (c1 = 0; c1 <= ub1; c1++) {
      {
        for (c2 = 0; c2 <= (_PB_NJ + -1); c2++) {
          {
            for (c5 = 0; c5 <= (_PB_NK + -1); c5++) {
              {
                tmp[c1][c2] += alpha * A[c1][c5] * B[c5][c2];
              }
            }
#pragma ivdep
#pragma vector always
#pragma simd
            for (c5 = 0; c5 <= (_PB_NL + -1); c5++) {
              {
                D[c1][c5] += tmp[c1][c2] * C[c2][c5];
              }
            }
          }
        }
      }
    }
  }
  if (((_PB_NJ >= 1) && (_PB_NK >= 1)) && (_PB_NL <= 0)) {
    ub1 = (_PB_NI + -1);
#pragma omp parallel for private(c5, c2) firstprivate(ub1)
    for (c1 = 0; c1 <= ub1; c1++) {
      {
        {
          for (c5 = 0; c5 <= (_PB_NK + -1); c5++) {
#pragma ivdep
#pragma vector always
#pragma simd
            for (c2 = 0; c2 <= (_PB_NJ + -1); c2++) {
              {
                tmp[c1][c2] += alpha * A[c1][c5] * B[c5][c2];
              }
            }
          }
        }
      }
    }
  }
  if (((_PB_NJ >= 1) && (_PB_NK <= 0)) && (_PB_NL >= 1)) {
    ub1 = (_PB_NI + -1);
#pragma omp parallel for private(c2, c5) firstprivate(ub1)
    for (c1 = 0; c1 <= ub1; c1++) {
      {
        for (c2 = 0; c2 <= (_PB_NJ + -1); c2++) {
          {
#pragma ivdep
#pragma vector always
#pragma simd
            for (c5 = 0; c5 <= (_PB_NL + -1); c5++) {
              {
                D[c1][c5] += tmp[c1][c2] * C[c2][c5];
              }
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
  int ni = NI;
  int nj = NJ;
  int nk = NK;
  int nl = NL;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(tmp,DATA_TYPE,NI,NJ,ni,nj);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NK,ni,nk);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NK,NJ,nk,nj);
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,NJ,NL,nj,nl);
  POLYBENCH_2D_ARRAY_DECL(D,DATA_TYPE,NI,NL,ni,nl);

  /* Initialize array(s). */
  init_array (ni, nj, nk, nl, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_2mm (ni, nj, nk, nl,
	      alpha, beta,
	      POLYBENCH_ARRAY(tmp),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nl,  POLYBENCH_ARRAY(D)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(tmp);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(D);

  return 0;
}
