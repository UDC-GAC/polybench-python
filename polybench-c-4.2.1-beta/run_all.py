import os,sys

selected_tests = [
"linear-algebra/blas/symm",
"linear-algebra/blas/syr2k",
"linear-algebra/blas/gemm",
"linear-algebra/blas/gesummv",
"linear-algebra/blas/trmm",
"linear-algebra/blas/syrk",
"linear-algebra/blas/gemver",
"linear-algebra/solvers/durbin",
"linear-algebra/solvers/ludcmp",
"linear-algebra/solvers/cholesky",
"linear-algebra/solvers/gramschmidt",
"linear-algebra/solvers/trisolv",
"linear-algebra/solvers/lu",
"linear-algebra/kernels/atax",
"linear-algebra/kernels/doitgen",
"linear-algebra/kernels/2mm",
"linear-algebra/kernels/mvt",
"linear-algebra/kernels/3mm",
"linear-algebra/kernels/bicg",
"stencils/seidel-2d",
"stencils/jacobi-2d",
"stencils/jacobi-1d",
"stencils/heat-3d",
"stencils/adi",
"stencils/fdtd-2d",
"datamining/covariance",
"datamining/correlation",
"medley/deriche",
"medley/floyd-warshall",
"medley/nussinov",
]

compiler = "gcc"
#compiler = "icc"
execution_count = 5
#opts = "-O1"
opts = "-O3 -march=skylake"
#opts = sys.argv[1]
#opts = "-O3 -march=skylake -fno-tree-vectorize -fno-tree-loop-vectorize"
#opts = "-O0"
cwd = os.getcwd()
for t in selected_tests:
    print("Running ", t.replace("/","."))
    print("Interpreter: C" )
    print("Options:")
    print("(%s,%d)"%("iterations",execution_count))
    print("(%s,%s)"%("'POLYBENCH_PAPI'",'True'))
    sys.stdout.flush()
    os.chdir( cwd + "/" + t )

    prog_name = t.split("/")[-1]
    src_base = prog_name + ".c"
    src_opt_pluto = prog_name + ".pocc.c"
    src_opt  = prog_name + "-prevector.pocc.c"

#    # Compile base
#    os.system( "%s %s -DPOLYBENCH_PAPI -I %s/utilities %s %s/utilities/polybench.c -lm -lpapi" % (compiler,opts,cwd,src_base,cwd) )
#
#    os.system( "cpupower frequency-set -d 3.7GHz" ) # Fix cpu frequency
#    os.system( "cpupower frequency-set -u 3.7GHz" ) # Fix cpu frequency
#
#    print( "Baseline:" )
#    sys.stdout.flush()
#    for i in range( execution_count ):
#        os.system( "./a.out" )
#
#    os.system( "cpupower frequency-set -d 800MHz" )
#    os.system( "cpupower frequency-set -u 4.7GHz" )  # Fix cpu frequency
#
#    continue
    # Compile opt
    ret = os.system( "~/workspace/pocc-devel/bin/pocc --pluto --pluto-prevector --vectorizer --pragmatizer --no-candl %s -o %s >/dev/null" % (src_base,src_opt) )
    if ret > 0:
        sys.exit(0)
    continue#
    os.system( "%s %s -DPOLYBENCH_PAPI -I %s/utilities %s %s/utilities/polybench.c -lm -lpapi" % (compiler,opts,cwd,src_opt,cwd) )

    os.system( "cpupower frequency-set -d 3.7GHz" ) # Fix cpu frequency
    os.system( "cpupower frequency-set -u 3.7GHz" ) # Fix cpu frequency

    print( "Pluto:" )
    sys.stdout.flush()
    for i in range( execution_count ):
        os.system( "./a.out" )

    os.system( "cpupower frequency-set -d 800MHz" )
    os.system( "cpupower frequency-set -u 4.7GHz" )  # Fix cpu frequency
