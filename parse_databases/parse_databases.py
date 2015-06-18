def read_decoy_database(filename):
    import csv
    sum_formulae = {}
    with open(filename) as filein_db:
        data_in = csv.reader(filein_db)
        for line in data_in:
            if line == []: #catch empty line(s)
                continue
            
            db_id= line[0]
            name= line[3]
            sf= line[2]
            mw = line[1]
            if float(mw) > 0: # useful filter for incomplete entries, not robust
                if sf in sum_formulae:
                    sum_formulae[sf]['name'].append(name)
                    sum_formulae[sf]['db_id'].append(db_id)
                else:
                    sum_formulae[sf] = {}
                    sum_formulae[sf]['name'] = [name]
                    sum_formulae[sf]['db_id'] = [db_id]
                    sum_formulae[sf]['mw'] = mw
        return sum_formulae

def read_kegg_compounds(filename):
    import csv
    sum_formulae = {}
    with open(filename) as filein_db:
        data_in = csv.reader(filein_db)
        for line in data_in:
            if line == []: #catch empty line(s)
                continue
            db_id = line[0]
            name= line[1]
            sf= line[2]
            mw= line[3]
            if float(mw) > 0: # useful filter for incomplete entries, not robust
                if sf in sum_formulae:
                    sum_formulae[sf]['name'].append(name)
                    sum_formulae[sf]['db_id'].append(db_id)
                else:
                    sum_formulae[sf] = {}
                    sum_formulae[sf]['name'] = [name]
                    sum_formulae[sf]['db_id'] = [db_id]
                    sum_formulae[sf]['mw'] = mw
        return sum_formulae
    
def read_hmdb_compounds(filename):
    import csv
    sum_formulae = {}
    with open(filename) as filein_db:
        data_in = csv.reader(filein_db,delimiter='\t')
        for line in data_in:
            if line == []: #catch empty line(s)
                continue
            db_id,mw,sf,name = line
            if mw=='':
                mw = '0'
            if float(mw) > 0: # useful filter for incomplete entries, not robust
                if sf in sum_formulae:
                    sum_formulae[sf]['name'].append(name)
                    sum_formulae[sf]['db_id'].append(db_id)
                else:
                    sum_formulae[sf] = {}
                    sum_formulae[sf]['name'] = [name]
                    sum_formulae[sf]['db_id'] = [db_id]
                    sum_formulae[sf]['mw'] = mw
    return sum_formulae

def read_pubchem_properties(filename):
    # download from pubchem - delete first blank line, add one blank line to end
    sum_formulae = {}
    l_count = 0
    with open(filename) as filein_db:
        for line in filein_db:
            if line == '\n': 
                if sf in sum_formulae:
                    sum_formulae[sf]['name'].append(name)
                    sum_formulae[sf]['db_id'].append(db_id)
                else:
                    sum_formulae[sf] = {}
                    sum_formulae[sf]['name'] = [name]
                    sum_formulae[sf]['db_id'] = [db_id]
                    sum_formulae[sf]['mw'] = mw
                l_count = 0
            else:
                if l_count==0: # name
                    name = line.split('.',1)[1]
                elif l_count==1: # formula is last string
                    info= line.strip().split() #remove whitespace, split on whitespace
                    try:
                        mw = info[info.index('MW:')+1]
                    except:
                        mw = ''
                    try:
                        sf = info[info.index('MF:')+1]
                        if '+' in sf or '-' in sf:
                            print 'bailing on charged molecule {} '.format(sf)
                            name=''
                            sf=''
                            mw=''
                            db_id=''
                            continue
                    except:
                        print 'failed for {}'.format(name)
                        continue
                elif l_count==2: # ID
                    db_id = line.strip().split()[1]
                l_count +=1
        sum_formulae.pop('',None)
    return sum_formulae
    