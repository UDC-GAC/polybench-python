import pandas as pd
import os, sys
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats

c_cols = ["cycles","inst","stalls","l1h","l1m","l2m","l3m","br","brmissp","scalard","p128d","p256d","loads","stores"]
abbrv_c_cols = ["cycles","inst","stalls","l1h","l1m","l2m","l3m","br","brmissp","loads","stores"]
py_cols= ["inst","cycles","stalls","l1m","l2m","l3m","br","brmissp","scalars","scalard","vecs","vecd","loads","stores"]

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
        print( "Usage: <whatever-fig.py> <results_folder> [output_image]" )
        print( "\tGenerates the requested file into [output_image], or opens it on the screen" )
        sys.exit(0)

    results_folder = sys.argv[1]
    if not os.path.isdir( results_folder ):
        print( "Error: %s is no directory" % (results_folder) )
        sys.exit(0)

    c_o0_path = results_folder + "/run_all_C_O0.out"
    c_o1_path = results_folder + "/run_all_C_O1.out"
    c_o3_path = results_folder + "/run_all_C_O3.out"
    c_o3_novec_path = results_folder + "/run_all_C_O3_novec.out"
    pypy_path = results_folder + "/pypy-run-all.out"
    map( nofile_error, [c_o0_path,c_o1_path,c_o3_path,c_o3_novec_path,pypy_path] )

    df_O0=parse(c_o0_path,cols=abbrv_c_cols)
    df_O1=parse(c_o1_path,cols=c_cols)
    df_O3=parse(c_o3_path,cols=c_cols)
    df_O3novec=parse(c_o3_novec_path,cols=c_cols)
    df_pypy=parse(pypy_path,cols=py_cols)

    df=pd.concat( [df_O0,df_O1,df_O3novec,df_O3,df_pypy], keys=["-O0","-O1","-O3 -fno-tree-vectorize","-O3","-O3","Python"], names=["gcc opts"] )
    df.reset_index().set_index( ["interpreter","benchmark","gcc opts","version"] ).sort_index(level=1).to_excel("/tmp/polybench_python.xls")

    o1_reuse = ['atax','bicg','mvt','gemver','gesummv','trisolv','deriche']
    on_reuse = ["correlation","covariance","2mm","3mm","gemm","symm","trmm","syrk","syr2k","cholesky","gramschmidt","lu","ludcmp","floyd-warshall","nussinov","durbin","doitgen"]
    stencils = ["fdtd-2d","jacobi-1d","jacobi-2d","heat-3d","seidel-2d","adi"]
    shape_1 = ["2mm","3mm","adi","cholesky","doitgen","lu","mvt","nussinov","seidel-2d","trisolv","trmm"]
    shape_2 = ["gesummv", "bicg", "syr2k", "ludcmp", "correlation", "covariance", "gramschmidt", "floyd-warshall", "symm"]
    shape_3 = ["jacobi-1d", "durbin", "atax", "gemver", "deriche", "syrk", "gemm", "fdtd-2d", "jacobi-2d", "heat-3d"]

    # Mask: remove PoCC and NumPy data points and the non-flattened version
    df=df.reorder_levels([1,2,0,3])
    mask = (df.index.get_level_values(3) == "baseline") | (df.index.get_level_values(3)=="list_flattened")
    df2= df[mask]
    df3 = df2.div(df2.xs( ["C","-O3 -fno-tree-vectorize","baseline"], level=[1,2,3] ))
    df4 = df3.stack().unstack([1,2,3])
    df4.columns = pd.Index( df4.columns )
    df4.rename( columns={ ('C', '-O1', 'baseline'): "gcc -O1", ('C', '-O3 -fno-tree-vectorize', 'baseline'): "gcc -O3 -fno-tree-vectorize", ('C', '-O3', 'baseline'):"gcc -O3", ('PyPy', 'Python', 'list'):"PyPy nested lists", ('PyPy', 'Python', 'list_flattened'): "PyPy flattened lists", ('PyPy','Python','load_elimination'): "PyPy with manual load elimination", ("C","-O0","baseline"): "gcc -O0" }, inplace=True )
    df4.drop( "gcc -O3 -fno-tree-vectorize", axis=1, inplace=True )
    df5 = df4.unstack().swaplevel(axis=1)

    sel_cols = ["inst","loads","stores","br"]
    df6=df5[sel_cols].apply( scipy.stats.gmean ).unstack()
#    df6['IPC of baseline']=(df2.inst/df2.cycles).loc[shape_1].xs("-O3 -fno-tree-vectorize",level=2).reset_index([1,2],drop=True)

    sns.set_context( "paper" )
    df6.plot(kind="bar",logy=False, colormap=sns.light_palette( "crimson", as_cmap=True ), edgecolor="k", rot=30 )
    plt.ylabel( "Normalized axis" )
    plt.axhline( 1.0, linestyle="--", color="k", linewidth=.5, label="gcc -O3 -fno-tree-vectorize" )
    plt.legend()
    if len(sys.argv) > 2:
        plt.savefig( sys.argv[2], bbox_inches='tight' )
    else:
        plt.show()

