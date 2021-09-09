def check_raw(motif,i):
    return [d for d in motif if d['position'] == i]

def analyze_index(Raw_Results):

    indexloss=[]
    indextomajor=[]
    indextominor=[]
    indextounique=[]
    for i in range (0,len(analyse.results)-1):
        t=check_raw(Raw_Results,i)

        if any(not d['hits'] for d in t):
            indexloss.append(i+1)
        if any('Major' in d['tag'] for d in t):
            indextomajor.append(i+1)
        if any('Minor' in d['tag'] for d in t):
            indextominor.append(i+1)
        if any('Unique' in d['tag'] for d in t):
            indextounique.append(i+1)

    return indexloss,indextomajor,indextominor,indextounique

def analyze_majors(Raw_Results):
    majorloss=[]
    majortoindex=[]
    majortominor=[]
    majortounique=[]

    for i in range (0,len(analyse.results)):
        t=check_raw(Raw_Results,i)

        if any(not d['hits'] for d in t):
            majorloss.append(i+1)
        if any('Index' in d['tag'] for d in t):
            majortoindex.append(i+1)
        if any('Minor' in d['tag'] for d in t):
            majortominor.append(i+1)
        if any('Unique' in d['tag'] for d in t):
            majortounique.append(i+1)

    return majorloss,majortoindex,majortominor,majortounique

def analyze_minors(Raw_Results):

    minorloss=[]
    minortoindex=[]
    minortomajor=[]
    minortounique=[]

    for i in range (0,len(analyse.results)):
        t=check_raw(Raw_Results,i)

        if any(not d['hits'] for d in t):
            minorloss.append(i+1)
        if any('Index' in d['tag'] for d in t):
            minortoindex.append(i+1)
        if any('Major' in d['tag'] for d in t):
            minortomajor.append(i+1)
        if any('Unique' in d['tag'] for d in t):
            minortounique.append(i+1)

    return minorloss,minortoindex,minortomajor,minortounique

def analyze_uniques(Raw_Results):

    uniqueloss=[]
    uniquetoindex=[]
    uniquetomajor=[]
    uniquetominor=[]
    for i in range (0,len(analyse.results)):
        t=check_raw(Raw_Results,i)

        if any(not d['hits'] for d in t):
            uniqueloss.append(i+1)
        if any('Index' in d['tag'] for d in t):
            uniquetoindex.append(i+1)
        if any('Major' in d['tag'] for d in t):
            uniquetomajor.append(i+1)
        if any('Minor' in d['tag'] for d in t):
            uniquetominor.append(i+1)
    return uniqueloss,uniquetoindex,uniquetomajor,uniquetominor

def gain(motif):
    gain=[]
    for x in range (1,len(analyse.results)):
        a = []
        b = []

        if analyse.results[x-1].support == 0 and analyse.results[x].support != 0:#no motif at the previous pos thus a gain
            for i in analyse.results[x].variants:
                if i.motif_long == motif:
                    gain.append(x+1)
            continue

        if analyse.results[x].support == 0:
            continue

        for i in analyse.results[x-1].variants: #get all the sequences in the previous pos, drop the first Seq
            a.append(i.sequence[1::])

        for i in analyse.results[x].variants: #get all the sequences related to motif, drop the last seq
            if i.motif_long == motif:
                b.append(i.sequence[0:-1])

        c=list(set(b) - set(a)) #get elements which are in b but not in a (previous pos)

        if c: #if there is any differense that means there is an alien seq in the motif
            gain.append(x+1)
    return gain

#Mergers, they either lose their motif and merge with the higher motif or
#there is a switch in a higher position, and the replacing pattern is coming from lower motif
#that indicates merging. be aware, a replacing pattern variant still can stay in the lower motif. Does not have to lose its motif.

#just like split.
#splits, while preserving its original motif, they either replace a motif or they spit a variant to the lower motif.

# while mergers increases the frequency, splits means a drop in the frequency.
#(unless there is a sudden dramatic increase or drop in sup)

def i_splits_mergers(Filtered_Indexes,
                    Raw_Indexes,
                    Raw_Majors,
                    Raw_Minors,
                    majortoindex,
                    minortoindex,
                    uniquetoindex):

    mtoi,mitoi,utoi=[],[],[]
    merge_mtoi,merge_mitoi,merge_utoi=[],[],[]

    for i in Filtered_Indexes:

        multi_in = [d for d in Raw_Indexes if d['position'] == i-1 and ' x ' in d['tag']] #get all motif changes related to pos//since we are only filtering only one at a position

        if len(multi_in) < len([d for d in Raw_Indexes if d['position'] == i-1]): #if there are 4 indexes, and 3 of them has switches, no need to look what became an index
            continue

        if i in majortoindex:
            mtoi.append(i)
            continue

        if check_raw(Raw_Majors,i-1):
            if any('Index' in str(d['hits'].keys()) for d in check_raw(Raw_Majors,i-1)):
                merge_mtoi.append(i)
                continue

        if i in minortoindex:
            mitoi.append(i)
            continue

        if check_raw(Raw_Minors,i-1):
            if any('Index' in str(d['hits'].keys()) for d in check_raw(Raw_Minors,i-1)):
                merge_mitoi.append(i)

        if i in uniquetoindex:
            utoi.append(i)
            continue

    for x in majortoindex:
        if x not in mtoi:
            merge_mtoi.append(x)
    for x in set(minortoindex):
        if x not in mitoi:
            merge_mitoi.append(x)
    for x in uniquetoindex:
        if x not in utoi:
            merge_utoi.append(x)
    return mtoi,mitoi,utoi,merge_mtoi,merge_mitoi,merge_utoi

def m_splits_mergers(Filtered_Majors,
                    Raw_Majors,
                    Raw_Indexes,
                    Raw_Minors,
                    indextomajor,
                    minortomajor,
                    uniquetomajor):

    split_itom,itom,mitom,merge_mitom,utom,merge_utom=([] for i in range(6))
    for i in Filtered_Majors:
        multi_in = [d for d in Raw_Majors if d['position'] == i-1 and ' x ' in d['tag']]

        if len(multi_in) < len([d for d in Raw_Majors if d['position'] == i-1]):
            continue

        if i in indextomajor:
            itom.append(i)
            continue
        if check_raw(Raw_Indexes,i-1):
            if any('Major' in str(d['hits'].keys()) for d in check_raw(Raw_Indexes,i-1)):
                split_itom.append(i)
                continue
        if i in minortomajor:
            mitom.append(i)
            continue
        if check_raw(Raw_Minors,i-1):
            if any('Major' in str(d['hits'].keys()) for d in check_raw(Raw_Minors,i-1)):
                merge_mitom.append(i)
                continue
        if i in uniquetomajor:
            utom.append(i)
            continue

        #No need below for unique cause unique has to come out
        #if check_raw('Unique',i-1):
            #if any('Major' in str(d['hits'].keys()) for d in check_raw('Unique',i-1)):
                #utom.append(i)
                #continue

    for x in minortomajor:
        if x not in mitom:
            merge_mitom.append(x)
    for x in uniquetomajor:
        if x not in utom:
            merge_utom.append(x)

    return split_itom,itom,mitom,merge_mitom,utom,merge_utom

def mi_splits_mergers(analyse,indextominor,majortominor,uniquetominor):

    split_itomi,split_mtomi,mtomi,utomi,merge_utomi,migain=([] for i in range(6))
    for x in range (1,len(analyse.results)):
        a = {'Index':[],'Major':[],'Minor':[],'Unique':[]}
        b = []

        if analyse.results[x-1].support == 0 and analyse.results[x].support != 0:#no motif at the previous pos thus a gain
            for i in analyse.results[x].variants:
                if i.motif_long == 'Minor':
                    migain.append(x+1)
            continue
        if analyse.results[x].support == 0:
            continue

        for i in analyse.results[x-1].variants: #get all the sequences in the previous pos, drop the first Seq and get the motif
            a[i.motif_long].append(i.sequence[1::])

        for i in analyse.results[x].variants: #get all the sequences related to motif, drop the last seq
            if i.motif_long == 'Minor':
                if i.sequence[0:-1] not in b:
                    b.append((i.sequence[0:-1]))

        a1=a['Minor']

        diff=list(set(b) - set(a1))

        if diff:
            for i in diff:
                if i in a['Index']:
                    if x+1 in indextominor:
                        continue
                    else:
                        split_itomi.append(x+1)
                elif i in a['Major']:
                    if x+1 in majortominor:
                        continue
                    else:
                        split_mtomi.append(x+1)
                elif i in a['Unique']:
                    #if x+1 in uniquetominor:
                    utomi.append(x+1)
                else:
                    migain.append(x+1)

    for x in uniquetominor:
        if x not in utomi:
            merge_utomi.append(x)
    utomi=list(set(utomi))
    merge_utomi=list(set(merge_utomi))
    migain=list(set(migain))

    return split_itomi,split_mtomi,mtomi,utomi,merge_utomi,migain

def u_splits_mergers(analyse,indextounique,majortounique,minortounique):

    split_itou,split_mtou,split_mitou,ugain=([] for i in range(4))

    for x in range (1,len(analyse.results)):
        a = {'Index':[],'Major':[],'Minor':[],'Unique':[]}
        b = []

        if analyse.results[x-1].support == 0 and analyse.results[x].support != 0:#no motif at the previous pos thus a gain
            for i in analyse.results[x].variants:
                if i.motif_long == 'Unique':
                    ugain.append(x+1)
            continue
        if analyse.results[x].support == 0:
            continue

        for i in analyse.results[x-1].variants: #get all the sequences in the previous pos, drop the first Seq and get the motif
            a[i.motif_long].append(i.sequence[1::])

        for i in analyse.results[x].variants: #get all the sequences related to motif, drop the last seq
            if i.motif_long == 'Unique':
                if i.sequence[0:-1] not in b:
                    b.append((i.sequence[0:-1]))

        a1=a['Unique']

        diff=list(set(b) - set(a1))

        if diff:
            for i in diff:
                if i in a['Index']:
                    if x+1 in indextounique:
                        continue
                    else:
                        split_itou.append(x+1)
                elif i in a['Major']:
                    if x+1 in majortounique:
                        continue
                    else:
                        split_mtou.append(x+1)
                elif i in a['Minor']:
                    if x+1 in minortounique:
                        continue
                    else:
                        split_mitou.append(x+1)
                else:
                    ugain.append(x+1)

    split_mitou=list(set(split_mitou))
    ugain=list(set(ugain))

    return split_itou,split_mtou,split_mitou,ugain
