import argparse
import SVDB_overlap_module
import SVDB_build_module
import SVDB_query_module
import SVDB_hist_module
import SVDB_purge_module
import SVDB_merge_vcf_module
import SVDB_query_module_SQLITE


if __name__ == '__main__':
    parser = argparse.ArgumentParser("""SVDB, use the build module to construct databases, use the query module to query the database usign vcf files, or use the hist module to generate histograms""",add_help=False)
    parser.add_argument('--build'       , help="create a db", required=False, action="store_true")
    parser.add_argument('--hist'        , help="generate histograms o the performance of a db", required=False, action="store_true")
    parser.add_argument('--query'       , help="query a db", required=False, action="store_true")
    parser.add_argument('--purge'       , help="remove entries from a database", required=False, action="store_true")
    parser.add_argument('--merge'       , help="merge similar structural variants within a vcf file", required=False, action="store_true")
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
        parser.add_argument('--ci', help="overides overlap and bnd_distance,determine hits based on the confidence interval of the position fo the variants(0 if no CIPOS or CIEND is vailable)", required=False, action="store_true")
        args= parser.parse_args()
        if args.folder or args.files:
            SVDB_hist_module.main(args)
        else:
            parser.print_help()
    elif args.query:
        parser = argparse.ArgumentParser("""SVDB: query module""")
        parser.add_argument('--query', help="query a db", required=False, action="store_true")
        parser.add_argument('--query_vcf', type=str, help="a vcf used to query the db", required = True)
        parser.add_argument('--db'        , type=str,  help="path to a SVDB db vcf ")
        parser.add_argument('--sqdb'        , type=str,  help="path to a SVDB sqlite db")
        parser.add_argument('--hit_tag'        , type=str,default="OCC",help="the tag used to describe the number of hits within the info field of the output vcf(default=OCC)")
        parser.add_argument('--frequency_tag'        , type=str,default="FRQ",help="the tag used to describe the frequency of the variant(defualt=FRQ)")
        parser.add_argument('--prefix', type=str,default=None ,help="the prefix of the output file, default = print to stdout")
        parser.add_argument('--bnd_distance', type=int,default= 10000,help="the maximum distance between two similar precise breakpoints(default = 10000)")
        parser.add_argument('--overlap', type=float, default = 0.6,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.6")
        parser.add_argument('--no_var',help="count overlaping variants of different type as hits in the db", required=False, action="store_true")
        parser.add_argument('--invert', help="invert the sorting order so that high frequency samples are present on top of the output vcf", required=False, action="store_true")
        parser.add_argument('--ci', help="overides overlap and bnd_distance,determine hits based on the confidence interval of the position fo the variants(0 if no CIPOS or CIEND is vailable)", required=False, action="store_true")
        args = parser.parse_args()

        if(args.db):
            SVDB_query_module.main(args)
        elif(args.sqdb):
            args.db=args.sqdb
            SVDB_query_module_SQLITE.main(args)
        else:
            print "invalid db option, choose --db to use the vcf db or sqdb to use the sqlite db"
    elif(args.build):
        parser = argparse.ArgumentParser("""SVDB: build module""")
        parser.add_argument('--build'       , help="create a db", required=False, action="store_true")
        parser.add_argument('--no_merge'       , help="skip the merging of variants", required=False, action="store_true")
        parser.add_argument('--ci', help="overides overlap and bnd_distance,determine hits based on the confidence interval of the position fo the variants(0 if no CIPOS or CIEND is vailable)", required=False, action="store_true")
        parser.add_argument('--bnd_distance', type=int,default= 2500,help="the maximum distance between two similar precise breakpoints(default = 2500)")
        parser.add_argument('--overlap', type=float, default = 0.8,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.8")
        parser.add_argument('--files'        , type=str, nargs='*', help="create a db using the specified vcf files(cannot be used with --folder)")
        parser.add_argument('--folder', type=str, help="create a db using all the vcf files in the folders")
        parser.add_argument('--prefix', type=str,default="SVDB" ,help="the prefix of the output file, default = SVDB")
        args = parser.parse_args()
        if (args.files and args.folder):
            print("ERROR: only one DB build input source may be selected");
            quit()
		#merging will be impossible
        if args.no_merge:
            args.overlap=float("inf")
            args.bnd_distance=-1
        if args.folder or args.files:
            SVDB_build_module.main(args)
        else:
            print("error, use files or folder to provide input for the database creation algorithm")

    elif args.purge:
        parser = argparse.ArgumentParser("""SVDB: purge module""")
        parser.add_argument('--purge', help="query a db", required=False, action="store_true")
        parser.add_argument('--vcf', type=str, help="remove all the entries within the DB that overlaps with any variant in the vcf")
        parser.add_argument('--samples', nargs='*',type=str, help="remove the given samples from the database")
        parser.add_argument('--file',type=str, help="remove the given samples from the database, each sample id is written, in a file, one line per sample id")
        parser.add_argument('--db'        , type=str, required=True, help="path to a SVDB vcf")
        parser.add_argument('--bnd_distance', type=int,default= 10000,help="the maximum distance between two similar precise breakpoints(default = 10000)")
        parser.add_argument('--overlap', type=float, default = 0.6,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.6")
        parser.add_argument('--ci', help="overides overlap and bnd_distance,determine hits based on the confidence interval of the position fo the variants(0 if no CIPOS or CIEND is vailable)", required=False, action="store_true")
        args= parser.parse_args()
        if args.samples or args.vcf or args.file:
            SVDB_purge_module.main(args)
        else:
            print("use the samples or vcf option to remove entries from the db")

    elif args.merge:
        parser = argparse.ArgumentParser("""SVDB: vcf_merge module""")
        parser.add_argument('--merge', help="merge structural variants", required=False, action="store_true")
        parser.add_argument('--vcf', nargs='*', type=str, help="remove all the entries within the DB that overlaps with any variant in the vcf",required=True)
        parser.add_argument('--bnd_distance', type=int,default= 2000,help="the maximum distance between two similar precise breakpoints(default = 2000)")
        parser.add_argument('--overlap', type=float, default = 0.95,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.95")
        parser.add_argument('--ci', help="overides overlap and bnd_distance,merge based on the confidence interval of the position fo the variants(0 if no CIPOS or CIEND is vailable)", required=False, action="store_true")
        parser.add_argument('--no_intra', help="no merging of variants within the same vcf", required=False, action="store_true")
        parser.add_argument('--no_var', help="variants of different type will be merged", required=False, action="store_true")
        parser.add_argument('--pass_only', help="merge only variants labeled PASS", required=False, action="store_true")
        args= parser.parse_args()        
        SVDB_merge_vcf_module.main(args)
    else:
        parser.print_help()

