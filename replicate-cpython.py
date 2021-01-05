#! /usr/bin/python


########################### Edit as needed
iterations = 5
dataset_size = "MINI" # This is the standard PolyBench size, used in the experiments in the paper
c_compiler = "/usr/bin/gcc"
O0_options = "-O0"
O1_options = "-O1"
O3_options = "-O3 -march=skylake" # Change -march as needed
O3_novec_options = O3_options + " -fno-tree-vectorize -fno-tree-loop-vectorize"
pypy_path = "pypy3"
cpython_path = "python3"
########################### The rest should take care of itself


import sys, os

try:
    import tqdm
    loop_iterator = tqdm.tqdm
    pipe_operator = ">>"
    stdout_file = open( "/dev/null", "w" )
    have_tqdm = True
    feedback_func = lambda t, x: t.set_description(x) and t.refresh()
except ModuleNotFoundError: 
    loop_iterator = list
    pipe_operator = "| tee -a"
    stdout_file = sys.stdout
    have_tqdm = False
    feedback_func = lambda t,x: print(x) and sys.stdout.flush()

all_tests = [
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

load_elimination_tests = [
"linear-algebra/blas/trmm",
"linear-algebra/solvers/cholesky",
"linear-algebra/solvers/trisolv",
"linear-algebra/solvers/lu",
"linear-algebra/kernels/doitgen",
"linear-algebra/kernels/2mm",
"linear-algebra/kernels/mvt",
"linear-algebra/kernels/3mm",
"stencils/adi",
"medley/nussinov",
]

def run_c( benchmark, results_path, version, compiler_options ):
    global benchmark_bar, c_compiler, root_folder, polybench_folder, c_polybench_options, iterations, pipe_operator, have_tqdm

    feedback_func( benchmark_bar, "Running %s/C %s %s" % (benchmark,c_compiler,compiler_options) )
    os.chdir( root_folder + polybench_folder + benchmark )
    c_file = benchmark.split("/")[-1] + ".c"
    os.system( "echo 'Running %s' >> %s" % (benchmark.replace("/","."),results_path) )
    os.system( "echo 'Interpreter: C' >> %s" % (results_path) ) 
    os.system( "echo 'Options:' >> %s" % (results_path) ) 
    os.system( "echo '(iterations,%d)' >> %s" % (iterations,results_path) ) 
    os.system( "echo '(\'POLYBENCH_PAPI\',True)' >> %s" % (results_path) ) 
    os.system( "echo 'Version=%s' >> %s" % (version,results_path) )
    os.system( "%s %s %s -I%s %s %s -lm -lpapi" % (c_compiler,compiler_options,c_polybench_options,root_folder + polybench_folder + "utilities", c_file, root_folder + polybench_folder + "utilities/polybench.c") )
    for i in range(iterations):
        os.system( "./a.out %s %s" %(pipe_operator, results_path) )
    os.remove( "a.out" )
    if have_tqdm: benchmark_bar.update(1)

def run_py( benchmark, results_path, version, py_interpreter, py_options ):
    global benchmark_bar, root_folder

    feedback_func( benchmark_bar, "Running %s/Python (%s) %s..." % (benchmark,version,py_interpreter) )
    os.chdir( root_folder )
    py_file = ("benchmarks/%s/%s.c" % (benchmark, benchmark.split("/")[-1])).replace("-","_")
    py_file = py_file.replace( "2mm", "_2mm" )
    py_file = py_file.replace( "3mm", "_3mm" )
    os.system( "echo 'Version=%s' >> %s" %(version,results_path) )
    os.system( "%s run-benchmark.py %s %s %s %s %s" % (py_interpreter,py_file,py_polybench_options,py_options,pipe_operator,results_path) )
    if have_tqdm: benchmark_bar.update(1)

c_polybench_options = "-DPOLYBENCH_PAPI -D%s_DATASET" % (dataset_size)
py_polybench_options = "--polybench-options POLYBENCH_PAPI --dataset-size %s --iterations %d" %(dataset_size,iterations)
polybench_folder = "polybench-c-4.2.1-beta/"

root_folder = os.getcwd() + "/"
results_path = root_folder + "results/"
regenerated_path = results_path + "regenerated/"

if not os.path.exists( results_path ):
    os.mkdir( results_path )

if not os.path.exists( regenerated_path ):
    os.mkdir( regenerated_path )

numpy_results_path = regenerated_path + "numpy-cpython-run-all.out"
try: os.remove( numpy_results_path )
except FileNotFoundError: pass

it = loop_iterator(all_tests)
for benchmark in it:
    feedback_func( it, "Regenerating results for %s ..." % (benchmark) )

    if have_tqdm:
        benchmark_bar = tqdm.trange(1, leave=False)
    else:
        benchmark_bar = None

    run_py( benchmark, numpy_results_path, "numpy", cpython_path, "--array-implementation 2" )

    if have_tqdm: benchmark_bar.close()
