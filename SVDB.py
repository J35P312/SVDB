import argparse
import SVDB_overlap_module
import SVDB_build_module
import SVDB_query_module
import SVDB_hist_module
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser("""SVDB, use the build module to construct databases, use the query module to query the database usign vcf files, or use the hist module to generate histograms""",add_help=False)
    parser.add_argument('--build'       , help="create a db", required=False, action="store_true")
    parser.add_argument('--hist'        , help="generate histograms o the performance of a db", required=False, action="store_true")
    parser.add_argument('--query'       , help="query a db", required=False, action="store_true")
    args, unknown = parser.parse_known_args()

    if args.hist:
        parser = argparse.ArgumentParser("""SVDB: query module""")
        parser.add_argument('--hist'        , help="generate histograms o the performance of a db", required=False, action="store_true")
        parser.add_argument('--similarity_matrix' , help ="output a csv table to stdout, the table indicates the nmber of similar svs of each input sample",required=False, action="store_true")
        parser.add_argument('--sample_hist', help ="accepts vcf files using folder or files. creates databases of size k for each given k m times. Reports an average frequency histogram per k as well as the number of unique samples per k,(database construction uses defaut db build settings)",required=False, action="store_true")
        parser.add_argument('--files'        , type=str, nargs='*', help="input vcf files(cannot be used with folder)")
        parser.add_argument('--folder', type=str, help="get all vcf files in the folder as input")
        parser.add_argument('--k', type=int,nargs='*',default=None, help="the number db sample sizes,default = n=10*i < samples")
        parser.add_argument('--n', type=int,default=100, help="the number of iterations,default=100(used with sample_hist)")
        parser.add_argument('--bnd_distance', type=int,default= 10000,help="the maximum distance between two similar precise breakpoints(default = 10000)")
        parser.add_argument('--overlap', type=float, default = 0.6,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.6")
        args= parser.parse_args()
        if args.folder or args.files:
            SVDB_hist_module.main(args)
        else:
            parser.print_help()
    elif args.query:
        parser = argparse.ArgumentParser("""SVDB: query module""")
        parser.add_argument('--query', help="query a db", required=False, action="store_true")
        parser.add_argument('--query_vcf', type=str, help="a vcf used to query the db")
        parser.add_argument('--db'        , type=str,  help="path to a SVDB vcf")
        parser.add_argument('--bnd_distance', type=int,default= 10000,help="the maximum distance between two similar precise breakpoints(default = 10000)")
        parser.add_argument('--overlap', type=float, default = 0.6,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.6")
        parser.add_argument('--no_var',help="count overlaping variants of different type as hits in the db", required=False, action="store_true")
        args = parser.parse_args()

        SVDB_query_module.main(args)
    elif(args.build):
        parser = argparse.ArgumentParser("""SVDB: build module""")
        parser.add_argument('--build'       , help="create a db", required=False, action="store_true")
        parser.add_argument('--bnd_distance', type=int,default= 5000,help="the maximum distance between two similar precise breakpoints(default = 5000)")
        parser.add_argument('--overlap', type=float, default = 0.8,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.8")
        parser.add_argument('--files'        , type=str, nargs='*', help="create a db using the specified vcf files(cannot be used with --folder)")
        parser.add_argument('--folder', type=str, help="create a db using all the vcf files in the folders")
        args = parser.parse_args()
        if (args.files and args.folder):
            print("ERROR: only one DB build input source may be selected");
            quit()
        SVDB_build_module.main(args)
    else:
        parser.print_help()

