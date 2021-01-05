import pandas as pd
import os, sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

c_cols = ["cycles","inst","stalls","l1h","l1m","l2m","l3m","br","brmissp","scalard","p128d","p256d","loads","stores"]
py_cols= ["inst","cycles","stalls","l1m","l2m","l3m","br","brmissp","scalars","scalard","vecs","vecd","loads","stores"]
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
        print( "Usage: <whatever-fig.py> <results_folder> [output_image]" )
        print( "\tGenerates the requested file into [output_image], or opens it on the screen" )
        sys.exit(0)

    results_folder = sys.argv[1]
    if not os.path.isdir( results_folder ):
        print( "Error: %s is no directory" % (results_folder) )
        sys.exit(0)

    c_o1_path = results_folder + "/run_all_C_O1.out"
    c_o3_path = results_folder + "/run_all_C_O3.out"
    c_o3_novec_path = results_folder + "/run_all_C_O3_novec.out"
    pypy_path = results_folder + "/pypy-run-all.out"
    map( nofile_error, [c_o1_path,c_o3_path,c_o3_novec_path,pypy_path] )

    df_O1=parse(c_o1_path,cols=c_cols)
    df_O3=parse(c_o3_path,cols=c_cols)
    df_O3novec=parse(c_o3_novec_path,cols=c_cols)
    try: df_pypy=parse(pypy_path,cols=py_cols)
    except ValueError: df_pypy=parse(pypy_path, cols=py_fp_cols2)

    df=pd.concat( [df_O1,df_O3novec,df_O3,df_pypy], keys=["-O1","-O3 -fno-tree-vectorize","-O3","Python"], names=["gcc opts"] )
    df.reset_index().set_index( ["interpreter","benchmark","gcc opts","version"] ).sort_index(level=1).to_excel("/tmp/polybench_python.xls")

    o1_reuse = ['atax','bicg','mvt','gemver','gesummv','trisolv','deriche']
    on_reuse = ["correlation","covariance","2mm","3mm","gemm","symm","trmm","syrk","syr2k","cholesky","gramschmidt","lu","ludcmp","floyd-warshall","nussinov","durbin","doitgen"]
    stencils = ["fdtd-2d","jacobi-1d","jacobi-2d","heat-3d","seidel-2d","adi"]
    shape_1 = ["2mm","3mm","adi","cholesky","doitgen","lu","mvt","nussinov","seidel-2d","trisolv","trmm"]
    shape_2 = ["gesummv", "bicg", "syr2k", "ludcmp", "correlation", "covariance", "gramschmidt", "floyd-warshall", "symm"]
    shape_3 = ["jacobi-1d", "durbin", "atax", "gemver", "deriche", "syrk", "gemm", "fdtd-2d", "jacobi-2d", "heat-3d"]
    sel_columns = ['cycles','inst','loads','l2m','l3m', 'stalls', 'integer']

    # Mask: remove PoCC and NumPy data points
    df = df[~(df.index.get_level_values(3).str.contains("pluto")) & (df.index.get_level_values(3) != "numpy")]
    df=df.reorder_levels([1,2,0,3])
    df2=pd.concat( [df.xs("list",level=3,drop_level=False),df.xs( "list_flattened", level=3, drop_level=False )] ).reset_index( level=[1,2], drop=True ).unstack().swaplevel(axis=1)

    groupby_dict = dict()
    for k in o1_reuse: groupby_dict[k] = "memory bound"
    for k in on_reuse: groupby_dict[k] = "compute bound"
    for k in stencils: groupby_dict[k] = "stencils"
    df3 = df2.groupby( groupby_dict ).sum()
    df3=df3.swaplevel(axis=1)
    int_column = df3.inst-df3.br-df3.scalard-df3.loads-df3.stores
    new_index = pd.MultiIndex( [["list","list_flattened"],["integer"]],  [[0,1],[0,0]] )
    int_column.columns = new_index
    df3=df3.swaplevel(axis=1)
    df3 = pd.concat( [df3, int_column], axis=1 )
    df4 = df3.list / df3.list_flattened

    mean_row = (df2.list / df2.list_flattened).mean()
    mean_row.name="aggregate"
    df5 = df4.append(mean_row)[sel_columns].dropna()

    sns.set_context( "paper" )
    df5[sel_columns].T.plot(kind="bar",logy=False, colormap=sns.light_palette( "crimson", as_cmap=True ), edgecolor="k", rot=0 )
    plt.axhline(1.0, linestyle="--", color="k", linewidth=.5 )
    plt.ylabel( "Normalized axis" )

    if len(sys.argv) > 2:
        plt.savefig( sys.argv[2], bbox_inches='tight' )
    else:
        plt.show()

