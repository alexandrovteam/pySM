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

def read_generic_csv(filename,idcol=0,namecol=1,mwcol=2,sfcol=3,header=1,verbose=False):
    import csv
    import re
    # some things lead to sf that cannot be parsed (mostly optional groups)
    # currently I just exclude these - a better database would provide combinations
    bad_groups = re.compile("(x)|(\,)|(\.n)|(\)n)|(R$)|(?:(R)[A-Z0-9])|(?:(X)[A-Z0-9])|(X$)")
    sum_formulae = {}
    with open(filename,'rU') as filein_db:
        data_in = csv.reader(filein_db,delimiter=',')
        # skip the headers
        for n in range(0,header):
            next(data_in, None)
        
        for line in data_in:
            if line == []: #catch empty line(s)
                continue
            db_id = line[idcol]
            sf = line[sfcol]
            sf=sf.replace(" ","")
            name = line[namecol]
            mw = line[mwcol]
            if mw=='':
                mw = '0'
            if '+' in sf or '-' in sf:
                print 'bailing on charged molecule {} '.format(sf)
                name=''
                sf=''
                mw=''
                db_id=''
                continue
            if any((sf =='',bad_groups.findall(sf)!=[])):
                if verbose:
                    print 'bad sf?: ID:{} name:{} sf:{}'.format(db_id,name,sf)
                continue
            if not sf in sum_formulae:
                sum_formulae[sf] = {'name':[], 'db_id':[], 'mw':[]}
            sum_formulae[sf]['name'].append(name)
            sum_formulae[sf]['db_id'].append(db_id)
            sum_formulae[sf]['mw'] = mw
    return sum_formulae

def read_helfrich_cyano_compounds(filename):
    import csv
    sum_formulae = {}
    with open(filename,'rU') as filein_db:
        data_in = csv.reader(filein_db,delimiter=',')

        for line in data_in:
            if line == []: #catch empty line(s)
                continue
            db_id,sf,name = line[0:3]
            if sf=='Formula': #skip header line
                continue
            mw=''
            if mw=='':
                mw = '0'
            if '+' in sf or '-' in sf:
                print 'bailing on charged molecule {} '.format(sf)
                name=''
                sf=''
                mw=''
                db_id=''
                continue
            if not sf in sum_formulae:
                sum_formulae[sf] = {'name':[], 'db_id':[], 'mw':[]}
            sum_formulae[sf]['name'].append(name)
            sum_formulae[sf]['db_id'] = db_id
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
    
