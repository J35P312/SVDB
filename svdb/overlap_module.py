from __future__ import absolute_import
#check the "overlap" of interchromosomaltranslocations
def precise_overlap(chrApos_query,chrBpos_query,chrApos_db,chrBpos_db,distance):    
    Adist=abs(chrApos_query-chrApos_db);
    Bdist=abs(chrBpos_query-chrBpos_db);
    if(Adist <= distance and Bdist <= distance ):
        return(True)

#use confidence intervals around the start and end position to determine similarity
def ci_overlap(chrApos_query,chrBpos_query,ciA_query,ciB_query,chrApos_db,chrBpos_db,ciA_db,ciB_db):

    query_A_upper = chrApos_query+abs(ciA_query[1])
    query_A_lower = chrApos_query-abs(ciA_query[0])

    query_B_upper = chrBpos_query+abs(ciB_query[1])
    query_B_lower = chrBpos_query-abs(ciB_query[0])

    #check if the position of the db sample is within the position of the query sample.
    if(query_A_upper >= chrApos_db  and chrApos_db >= query_A_lower and query_B_upper >= chrBpos_db and chrBpos_db >= query_B_lower):
        return(True) 
    db_A_upper = chrApos_db+abs(ciA_db[1])
    db_A_lower = chrApos_db-abs(ciA_db[0])

    db_B_upper = chrBpos_db+abs(ciB_db[1])
    db_B_lower = chrBpos_db-abs(ciB_db[0])

    if(db_A_upper >= chrApos_query  and chrApos_query >= db_A_lower and db_B_upper >= chrBpos_query and chrBpos_query >= db_B_lower):
        return(True)

#use confidence intervals around the start and end position to determine similarity, both confidence intervals must overlap both breakpoints
def ci_overlap_two_sided(chrApos_query,chrBpos_query,ciA_query,ciB_query,chrApos_db,chrBpos_db,ciA_db,ciB_db):

    query_A_upper = chrApos_query+abs(ciA_query[1])
    query_A_lower = chrApos_query-abs(ciA_query[0])

    query_B_upper = chrBpos_query+abs(ciB_query[1])
    query_B_lower = chrBpos_query-abs(ciB_query[0])

    #check if the position of the db sample is within the position of the query sample.
    first =(query_A_upper >= chrApos_db  and chrApos_db >= query_A_lower and query_B_upper >= chrBpos_db and chrBpos_db >= query_B_lower)

    db_A_upper = chrApos_db+abs(ciA_db[1])
    db_A_lower = chrApos_db-abs(ciA_db[0])

    db_B_upper = chrBpos_db+abs(ciB_db[1])
    db_B_lower = chrBpos_db-abs(ciB_db[0])

    second= (db_A_upper >= chrApos_query  and chrApos_query >= db_A_lower and db_B_upper >= chrBpos_query and chrBpos_query >= db_B_lower)
    if first and second:
        return(True)


#check if intrachromosomal vaiants overlap
def isSameVariation(chrApos_query,chrBpos_query,chrApos_db,chrBpos_db,ratio,distance): #event is in the DB, variation is the new variation I want to insert
    if abs(chrApos_query-chrApos_db) <= distance and abs(chrBpos_query-chrBpos_db) <= distance:
        if(chrApos_query < chrApos_db):
            region_start=chrApos_query;
            overlap_start=chrApos_db;
        else:
            region_start=chrApos_db;
            overlap_start=chrApos_query;

        if(chrBpos_query > chrBpos_db):
            region_end=chrBpos_query
            overlap_end=chrBpos_db

        else:
            region_end=chrBpos_db
            overlap_end=chrBpos_query
        try:
            event_ratio=float(overlap_end-overlap_start+1)/float(region_end-region_start+1)
        except:
            event_ratio=0;
        if event_ratio >= ratio:
            return(True)
        return None
    else:
        return None
def variant_overlap(chrA,chrB,chrApos_query,chrBpos_query,chrApos_db,chrBpos_db,ratio,distance):

    overlap = False
    if chrA == chrB:
        overlap= isSameVariation(chrApos_query,chrBpos_query,chrApos_db,chrBpos_db,ratio,distance)
    else:
        overlap = precise_overlap(chrApos_query,chrBpos_query,chrApos_db,chrBpos_db,distance)
    return(overlap)
