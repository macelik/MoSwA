def analyze_index(which,indexes,majors,minors,uniques):
    indexlost=[]
    majortoindex=[]
    minortoindex=[]
    uniquetoindex=[]
    indexgain=[]
    for i in indexes:# works v3 fixed
            if indexes[i]['hits'] == {}:
                    indexlost.append(i)
                    continue
            if i in majors:              
                if 'Index' in majors[i]['tag']:                    
                    majortoindex.append(i)
                    continue
            if i in minors:                
                if 'Index' in minors[i]['tag']:
                    minortoindex.append(i)
                    continue    
            if i in uniques:                
                if 'Index' in uniques[i]['tag']:                
                    uniquetoindex.append(i)
                    continue    
                else:
                    indexgain.append(i)
                    continue
            else:
                indexgain.append(i)

    return indexlost,majortoindex,minortoindex,uniquetoindex,indexgain

def average_pos(myd):
    t_count=0
    e_pos=[]
    for i in range(len(myd)):
        if myd[i]['sequences'] == []:
            e_pos.append(myd[i]['position'])
        for x in range(0,len(myd[i]['sequences'])):
            t_count+=myd[i]['sequences'][x]['count']

    average=t_count/len(myd)
    average="%.2f" % average
    return (average,e_pos)

def smthng():
    below=[]
    above=[]
    average_pos()
    for i in somethingelse:
        if indexes[i]['count'] > thold:
            above.append(i)
        else:
            below.append(i)
    return below,above

def analyze_major(which,indexes,majors,minors,uniques):
    majorlost,indextomajor,minortomajor,uniquetomajor,majorgain=[],[],[],[],[]
    for i in majors:
        if majors[i]['hits'] == {}:
            majorlost.append(i)
            continue
        if i in indexes:              
            if 'Major' in indexes[i]['tag']:                    
                indextomajor.append(i)
                continue
        if i in minors:                
            if 'Major' in minors[i]['tag']:
                minortomajor.append(i)
                continue    
        if i in uniques:                
            if 'Major' in uniques[i]['tag']:                
                uniquetomajor.append(i)
                continue    
            else:
                majorgain.append(i)
                continue
        else:
            majorgain.append(i)
    return majorlost,indextomajor,minortomajor,uniquetomajor,majorgain

def analyze_minor(which,indexes,majors,minors,uniques):
    minorlost,indextominor,majortominor,uniquetominor,minorgain=[],[],[],[],[]
    for i in minors:# works v3 fixed
        if minors[i]['hits'] == {}:
                minorlost.append(i)
                continue
        if i in indexes:              
            if 'Minor' in indexes[i]['tag']:                    
                indextominor.append(i)
                continue
        if i in majors:                
            if 'Minor' in majors[i]['tag']:
                majortominor.append(i)
                continue    
        if i in uniques:                
            if 'Minor' in uniques[i]['tag']:                
                uniquetominor.append(i)
                continue    
            else:
                minorgain.append(i)
                continue
        else:
            minorgain.append(i)
    return minorlost,indextominor,majortominor,uniquetominor,minorgain

def analyze_unique(which,indexes,majors,minors,uniques):
    uniquelost,indextounique,majortounique,minortounique,uniquegain=[],[],[],[],[]
    for i in uniques:# works v3 fixed
        if uniques[i]['hits'] == {}:
                uniquelost.append(i)
                continue
        if i in indexes:              
            if 'Unique' in indexes[i]['tag']:                    
                indextounique.append(i)
                continue
        if i in majors:                
            if 'Unique' in majors[i]['tag']:
                majortounique.append(i)
                continue    
        if i in minors:                
            if 'Unique' in minors[i]['tag']:                
                minortounique.append(i)
                continue    
            else:
                uniquegain.append(i)
                continue
        else:
            uniquegain.append(i)
    return uniquelost,indextounique,majortounique,minortounique,uniquegain