import pandas as pd
import sys, os
import matplotlib.pyplot as plt
import seaborn as sns

c_cols = ["cycles","inst","stalls","l1h","l1m","l2m","l3m","br","brmissp","scalard","p128d","p256d","loads","stores"]
abbrv_c_cols = ["cycles","inst","stalls","l1h","l1m","l2m","l3m","br","brmissp","loads","stores"]
py_cols= ["inst","cycles","stalls","l1m","l2m","l3m","br","brmissp","scalars","scalard","vecs","vecd","loads","stores"]
py_fp_cols= ["inst","cycles","stalls","l1m","l2m","l3m","br","brmissp","scalars_trash","scalard_trash","scalard","p128d","p256d","loads","stores"]
py_fp_cols2= ["inst","cycles","stalls","l1m","l2m","l3m","br","brmissp","scalard","p128d","p256d","loads","stores"]

def parse( path, cols=c_cols ):
    f = open( path, "r" )

    records = []
    options = {}
    for l in f:
        l = l.strip()
        if l.startswith( "Setting" ): continue
        if l.startswith( "Running" ):
            benchmark = l.split(".")[-1]
            parse_options = False
            options = {}
            iterations = 0
            continue
        if l.startswith( "Start time" ): continue
        if l.startswith( "Version=" ):
            version = l.split("=")[-1]
            if "iterations" in options: iterations = options['iterations']
            continue
        if l.startswith( "Interpreter:" ):
            interpreter = l.split( " " )[-1]
            continue
        if l.startswith( "Options:" ):
            parse_options = True
            continue
        if parse_options and l.startswith( "(" ):
            l = l.split("(")[-1].split(")")[0]
            param, value = l.split(",")
            if param[0] == "'": param = param[1:-1]
            try: options[param.strip()] = int(value.strip())
            except: options[param.strip()] = value.strip()
            if param == "iterations":
                iterations = int(value.strip())
            continue
        if parse_options and iterations > 0:
            if ('POLYBENCH_TIME' in options) and options['POLYBENCH_TIME'] == "True":
                measurements = [ float( l ) ]
            elif ('POLYBENCH_PAPI' in options) and options['POLYBENCH_PAPI'] == "True":
                measurements = list( map( float, l.split( " " ) ) )
            else: raise ValueError
            iterations -= 1
            records.append( [benchmark, interpreter, version, options] + measurements )
            continue
        raise ValueError

    df=pd.DataFrame(records,columns=["benchmark","interpreter","version","options"] + cols)
    df.benchmark=df.benchmark.str.replace("_2mm","2mm")
    df.benchmark=df.benchmark.str.replace("_3mm","3mm")
    df.benchmark=df.benchmark.str.replace("fdtd_2d","fdtd-2d")
    df.benchmark=df.benchmark.str.replace("heat_3d","heat-3d")
    df.benchmark=df.benchmark.str.replace("jacobi_1d","jacobi-1d")
    df.benchmark=df.benchmark.str.replace("jacobi_2d","jacobi-2d")
    df.benchmark=df.benchmark.str.replace("seidel_2d","seidel-2d")
    df.benchmark=df.benchmark.str.replace("floyd_warshall","floyd-warshall")
    df.drop( "options", axis=1, inplace=True )
#    gb = df.groupby( ["benchmark","interpreter","version"] )
#    df = gb.apply( lambda x: x.sort_values(by="cycles").iloc[0] )
#    df = df.groupby(level=range(5)).apply( lambda x: x.sort_values(by="cycles").iloc[0] ) # For each set of repetitions, retain the row with the best execution time
#    df = df[((gb.apply( lambda x: x-x.mean() ) / gb.transform("std")).abs().fillna(0) < 3).all(axis=1)]
#    df = df[(gb.apply( lambda x: x-x.mean() ).divide( gb.std() ).abs().fillna(0) <= 3).all( axis=1 ).values] # Discard observations which are more than 3*std away from the mean

    return df.groupby( ["benchmark","interpreter","version"] ).mean()
#    return df.reset_index(drop=True).set_index( ["benchmark","interpreter","version"] )

def nofile_error( path ):
    if not os.path.exists( path ):
        print( "Error: results file \"%s\" does not exist." % (path) )
        sys.exit( 0 )

if __name__ == "__main__":

    if len( sys.argv ) < 2:
        print( "Usage: <whatever-fig.py> <results_folder>" )
        print( "\tGenerates the appropriate image and tables into <results_folder>" )
        sys.exit(0)

    results_folder = sys.argv[1]
    if not os.path.isdir( results_folder ):
        print( "Error: %s is no directory" % (results_folder) )
        sys.exit(0)
    
    c_o3_path = results_folder + "/run_all_C_O3.out"
    c_o3_prevector_path = results_folder + "/run_all_C_O3_prevector.out"
    c_o3_novec_path = results_folder + "/run_all_C_O3_novec.out"
    pypy_path = results_folder + "/pypy-run-all.out"
    numpy_path = results_folder + "/numpy-cpython-run-all.out"
    map( nofile_error, [c_o3_path,c_o3_prevector_path,c_o3_novec_path,pypy_path,numpy_path] )

    df_O3=parse(c_o3_path,cols=c_cols)
    df_O3_prevector=parse(c_o3_prevector_path,cols=c_cols)
    df_O3novec=parse(c_o3_novec_path,cols=c_cols)
    try: df_pypy=parse(pypy_path,cols=py_cols)
    except ValueError: df_pypy = parse( pypy_path, cols=py_fp_cols2 )
    df_numpy=parse(numpy_path,cols=py_fp_cols2)

    df=pd.concat( [df_O3novec,df_O3,df_O3_prevector,df_pypy,df_numpy], keys=["-O3 -fno-tree-vectorize","-O3","-O3","Python","Python"], names=["gcc opts"] )
    df.reset_index().set_index( ["interpreter","benchmark","gcc opts","version"] ).sort_index(level=1).to_excel("/tmp/polybench_python.xls")

    o1_reuse = ['atax','bicg','mvt','gemver','gesummv','trisolv','deriche']
    on_reuse = ["correlation","covariance","2mm","3mm","gemm","symm","trmm","syrk","syr2k","cholesky","gramschmidt","lu","ludcmp","floyd-warshall","nussinov","durbin","doitgen"]
    stencils = ["fdtd-2d","jacobi-1d","jacobi-2d","heat-3d","seidel-2d","adi"]
    shape_1 = ["2mm","3mm","adi","cholesky","doitgen","lu","mvt","nussinov","seidel-2d","trisolv","trmm"]
    shape_2 = ["jacobi-1d", "durbin", "atax", "gemver", "deriche", "syrk", "gemm", "fdtd-2d", "jacobi-2d", "heat-3d"]
    shape_3 = ["gesummv", "bicg", "syr2k", "ludcmp", "correlation", "covariance", "gramschmidt", "floyd-warshall", "symm"]

    # Mask: remove PoCC and NumPy data points and the non-flattened version
    df=df.reorder_levels([1,2,0,3])

    # Fix jacobi-1d and durbin: list is list_flattened
    df.reset_index(level=3,inplace=True)
    for k in ['jacobi-1d','durbin']:
        row = df.loc[k,'PyPy']
        row[row.version=="list_pluto"].version = "list_flattened_pluto"
        row[row.version=="list"].version = "list_flattened"
    df.set_index("version",append=True,inplace=True)
    
    df2 = df.loc[['gemm','gramschmidt','syr2k',"seidel-2d","jacobi-2d"]]
    # Remove redundancy between numpy and numpy_wf. Keep best execution time.
#    df3 = df2.groupby( level=[0,1,2,3] ).apply( lambda x: x.sort_values(by="cycles").iloc[0] )
    df3 = df2.groupby( level=[0,1,2,3] ).apply( lambda x: x.sort_values(by="cycles").iloc[0] )
    df4 = df3/df3.xs(["C","-O3","baseline"], level=[1,2,3] )
    df5 = df4.cycles.unstack([1,2,3])
    df5.columns = pd.Index( df5.columns )
    df5.rename( columns={ ('C', '-O1', 'baseline'): "gcc -O1", ('C', '-O3 -fno-tree-vectorize', 'baseline'): "gcc -O3 -fno-tree-vectorize", ('C', '-O3', 'baseline'):"gcc -O3", ('PyPy', 'Python', 'list'):"PyPy nested lists", ('PyPy', 'Python', 'list_flattened'): "PyPy", ('PyPy','Python','load_elimination'): "PyPy with manual load elimination", ('CPython','Python','numpy'):'CPython/NumPy', ("C","-O3 -fno-tree-vectorize","pluto"): "PoCC --pluto / gcc -O3 -fno-tree-vectorize", ("C","-O3","pluto"): "pocc --pluto / gcc -O3", ("PyPy","Python","list_flattened_pluto"):"pocc --pluto / PyPy" }, inplace=True )
#    sel_rows = ["cycles", "inst", "l3m", "br", "loads","stores", "p256d"]
    sel_cols = ["pocc --pluto / gcc -O3", "PyPy", "pocc --pluto / PyPy", "CPython/NumPy"]


    sns.set_context( "paper" )
    df5[sel_cols].sort_index().plot(kind="bar",logy=True, colormap=sns.light_palette( "crimson", as_cmap=True), edgecolor="k", rot=0 )
    plt.ylabel( "Normalized execution cycles" )
    plt.axhline( 1.0, linestyle="--", color="k", linewidth=.5, label="gcc -O3" )
    plt.legend()
    plt.savefig( results_folder+"/fig1.pdf", bbox_inches='tight' )

