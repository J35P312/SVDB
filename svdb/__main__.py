import argparse
from . import overlap_module
from . import build_module
from . import query_module
from . import merge_vcf_module
from . import export_module

def main():
    version = "2.0.0"
    parser = argparse.ArgumentParser("""SVDB-{}, use the build module to construct databases, use the query module to query the database usign vcf files, or use the hist module to generate histograms""".format(version),add_help=False)
    parser.add_argument('--build'       , help="create a db", required=False, action="store_true")
    parser.add_argument('--query'       , help="query a db", required=False, action="store_true")
    parser.add_argument('--merge'       , help="merge similar structural variants within a vcf file", required=False, action="store_true")
    parser.add_argument('--export'       , help="export a database", required=False, action="store_true")
    args, unknown = parser.parse_known_args()

    if args.query:
        parser = argparse.ArgumentParser("""SVDB.{}: query module""".format(version))
        parser.add_argument('--query', help="query a db", required=False, action="store_true")
        parser.add_argument('--query_vcf', type=str, help="a vcf used to query the db", required = True)
        parser.add_argument('--db'        , type=str,  help="path to a SVDB db vcf ")
        parser.add_argument('--sqdb'        , type=str,  help="path to a SVDB sqlite db")
        parser.add_argument('--bedpedb'        , type=str,  help="path to a SV database of the following format chrA-posA-chrB-posB-type-count-frequency")
        parser.add_argument('--in_occ'        , type=str,help="The allele count tag, if used, this tag must be present in the INFO column of the input DB(usually set to AC or OCC)")
        parser.add_argument('--in_frq'        , type=str,help="The frequency count tag, if used, this tag must be present in the INFO column of the input DB(usually set to AF or FRQ)")
        parser.add_argument('--out_occ'        , type=str,default="OCC",help="the allle count tag, as annotated by SVDBvariant(defualt=OCC)")
        parser.add_argument('--out_frq'        , type=str,default="FRQ",help="the tag used to describe the frequency of the variant(defualt=FRQ)")
        parser.add_argument('--prefix', type=str,default=None ,help="the prefix of the output file, default = print to stdout")
        parser.add_argument('--bnd_distance', type=int,default= 10000,help="the maximum distance between two similar breakpoints(default = 10000)")
        parser.add_argument('--overlap', type=float, default = 0.6,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.6")
        parser.add_argument('--memory'       , help="load the database into memory: increases the memory requirements, but lowers the time consumption(may only be used with sqdb)", required=False, action="store_true")
        parser.add_argument('--no_var',help="count overlaping variants of different type as hits in the db", required=False, action="store_true")
        args = parser.parse_args()
        args.version=version

        if(args.db or args.sqdb or args.bedpedb):
            query_module.main(args)
        else:
            print ("invalid db option, choose --db to use the vcf db or sqdb to use the sqlite db")

    elif(args.build):
        parser = argparse.ArgumentParser("""SVDB-{}: build module""".format(version))
        parser.add_argument('--build'       , help="create a db", required=False, action="store_true")
        parser.add_argument('--files'        , type=str, nargs='*', help="create a db using the specified vcf files(cannot be used with --folder)")
        parser.add_argument('--folder', type=str, help="create a db using all the vcf files in the folders")
        parser.add_argument('--prefix', type=str,default="SVDB" ,help="the prefix of the output file, default = SVDB")
        args = parser.parse_args()
        args.version=version
        if (args.files and args.folder):
            print("ERROR: only one DB build input source may be selected");
            quit()

        if args.folder or args.files:
            build_module.main(args)
        else:
            print("error, use files or folder to provide input for the database creation algorithm")
    elif args.export:
        parser = argparse.ArgumentParser("""SVDB-{}: export module; export the variants of the SVDB sqlite database into a vcf file""".format(version))
        parser.add_argument('--export'       , help="create a db", required=False, action="store_true")
        parser.add_argument('--db'        , type=str, required=True, help="The SQLite database")
        parser.add_argument('--no_merge'       , help="skip the merging of variants, print all variants in the db to a vcf file", required=False, action="store_true")
        parser.add_argument('--bnd_distance', type=int,default= 2500,help="the maximum distance between two similar precise breakpoints(default = 2500)")
        parser.add_argument('--overlap', type=float, default = 0.8,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.8")
        parser.add_argument('--DBSCAN'       , help="use dbscan to cluster the variants", required=False, action="store_true")
        parser.add_argument('--epsilon'       ,type=float, default = 500, help="used together with --DBSCAN; sets the epsilon paramter(default = 500)", required=False)
        parser.add_argument('--min_pts'       ,type=int, default = 2, help="used together with 1--DBSCAN; sets the min_pts parameter(default = 2)", required=False)              
        parser.add_argument('--prefix', type=str,default="SVDB" ,help="the prefix of the output file, default = same as input")
        parser.add_argument('--memory'       , help="load the database into memory: increases the memory requirements, but lowers the time consumption", required=False, action="store_true")
        args = parser.parse_args()
        
        #merging will be impossible
        if args.no_merge:
            args.overlap=float("inf")
            args.bnd_distance=-1
        
        args.version=version        
        export_module.main(args)

    elif args.merge:
        parser = argparse.ArgumentParser("""SVDB-{}: vcf_merge module""".format(version))
        parser.add_argument('--merge', help="merge structural variants", required=False, action="store_true")
        parser.add_argument('--notag', help="Do not add the the VARID and set entries to the info field", required=False, action="store_true")
        parser.add_argument('--vcf', nargs='*', type=str, help="input vcf files, all input vcf files will be merged into one. Use the --prioriy flag to prioritize the callers/vcf files",required=True)
        parser.add_argument('--priority', type=str, help="prioritise the input files, using the following format --vcf caller1.vcf:2 caller2.vcf:1 --priority: 1,2")
        parser.add_argument('--bnd_distance', type=int,default= 2000,help="the maximum distance between two similar precise breakpoints(default = 2000)")
        parser.add_argument('--overlap', type=float, default = 0.95,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.95")
        parser.add_argument('--no_intra', help="no merging of variants within the same vcf", required=False, action="store_true")
        parser.add_argument('--no_var', help="variants of different type will be merged", required=False, action="store_true")
        parser.add_argument('--pass_only', help="merge only variants labeled PASS", required=False, action="store_true")
        parser.add_argument('--same_order', help="Across all input vcf files, the order of the sample columns are the same", required=False, action="store_true")
        args= parser.parse_args() 
        args.version=version       
        merge_vcf_module.main(args)
        
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
