#check the "overlap" of interchromosomaltranslocations
def precise_overlap(chrApos_query,chrBpos_query,chrApos_db,chrBpos_db,distance):    
    Adist=abs(chrApos_query-chrApos_db);
    Bdist=abs(chrBpos_query-chrBpos_db);
    if(Adist <= distance and Bdist <= distance ):
        return(True)

#check if intrachromosoma vaiants overlap
def isSameVariation(chrApos_query,chrBpos_query,chrApos_db,chrBpos_db,ratio): #event is in the DB, variation is the new variation I want to insert

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
