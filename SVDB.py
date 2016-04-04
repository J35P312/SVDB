import argparse
import SVDB_overlap_module
import SVDB_build_module
import SVDB_query_module
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser("""
    This scripts wants as input a databse containing varaition and counts [chrA chrB startA endA startB endB occurences]
    and a variation file produced by FindTranslocations or CNVnator. The script will output the variation file sorting the variation
    by number of occurences in the DB.
    """)
    parser.add_argument('--query_vcf', type=str, help="a vcf used to query the db")
    parser.add_argument('--db'        , type=str,  help="path to a SVDB vcf")
    parser.add_argument('--bnd_distance', type=int,default= 10000,help="the maximum distance between two similar precise breakpoints(default = 10000)")
    parser.add_argument('--overlap', type=float, default = 0.6,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.6")
    parser.add_argument('--no_var',help="count overlaping variants of different type as hits in the db", required=False, action="store_true")
    parser.add_argument('--build', type=str, help="create a db using all the vcf files in the folders")
    parser.add_argument('--files'        , type=str, nargs='*', help="creata a db using the specified vcf files(cannot be used with --build)")
    args = parser.parse_args()

    if (args.query_vcf or args.db):
        SVDB_query_module.main(args)
    elif(args.build):
        if (args.files and args.db):
            print("ERROR: only one DB build input source may be selected");
            quit()
        SVDB_build_module.main(args)
    else:
        print("no option selected, use --help to get more info")

