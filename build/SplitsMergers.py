import re
def view_json():
    vis={
        "MoSwa_Output":
        {
            "Report":{
                "kmer_length":"Null",
                "tresh_hold":"Null",
                "alignment_length":"Null",
                "highest_support": "Null",
                "average_support": "Null",
                "Positions_no_support":"Null",
                "Positions_low_support":"Null",
                "Total_switches":"Null",
                "Unique_switches":"Null",
                "Results":{
                    "Index_Motifs":{
                    "Index_Raw":{},
                    "Index_switches":{"ToMajor":{},
                                     "ToMinor":{},
                                     "ToUnique":{},
                                     "Lost":{},
                                     "Gain":{}},
                    "Index_Split":{"ToMajor":{},
                                  "ToMinor":{},
                                  "ToUnique":{}}
                    },
                    "Major_Motifs":{
                    "Major_Raw":{},
                    "Major_switches":{"ToIndex":{},
                                     "ToMinor":{},
                                     "ToUnique":{},
                                     "Lost":{},
                                     "Gain":{}},
                    "Major_Split":{"ToMinor":{},
                                  "ToUnique":{}},
                    "Major_Merger":{"ToIndex":{}}
                    },
                    "Minor_Motifs":{
                    "Minor_Raw":{},
                    "Minor_switches":{"ToIndex":{},
                                     "ToMajor":{},
                                     "ToUnique":{},
                                     "Lost":{},
                                     "Gain":{}},
                    "Minor_Split":{"ToUnique":{}},
                    "Minor_Merger":{"ToIndex":{},
                                   "ToMajor":{}}
                    },
                    "Unique_Motifs":{
                    "Unique_Raw":{},
                    "Unique_switches":{"ToIndex":{},
                                     "ToMajor":{},
                                     "ToMinor":{},
                                     "Lost":{},
                                     "Gain":{}},
                    "Unique_Merger":{"ToIndex":{},
                                   "ToMajor":{},
                                    "ToMinor":{}}
                    }

                },
                "PairWise_Alignments":{}
            }
        }
    }
    return vis

def insert_data(a):
    up=[]
    for x in a:
        d=dict()
        d['position']=x['position']
        d['seq']=x['seq']
        d['frequency']=x['count']
        d['tag']=x['tag']
        d['hits']=[]
        hits=[]
        count=0
        upd={}
        for i,j in x['hits'].items():
            count = count+1
            if count % 2 != 0:

                seq=i
                freq=j.split(" ")[-1]
                upd['position']=x['position']+1
                upd['seq']=seq
                upd['frequency']=freq
            else:
                tag=i
                upd['tag']=tag
                hits.append(upd)
                upd={}
        d['hits']=hits
        up.append(d)
    return up

def check_raw(motif,i):
    return [d for d in motif if d['position'] == i]

def analyze_index(Raw_Results,vis):

    indexloss=[]
    indextomajor=[]
    indextominor=[]
    indextounique=[]
    for i in range (1,len(analyse.results)):
        t=check_raw(Raw_Results,i)
        up=insert_data(t)
        vis["MoSwa_Output"]["Report"]["Results"]['Index_Motifs']['Index_Raw'][i]=up

        if any(not d['hits'] for d in t):
            indexloss.append(i+1)
            z=[d for i,d in enumerate(t) if not d['hits']]
            up=insert_data(z)
            vis["MoSwa_Output"]["Report"]["Results"]['Index_Motifs']['Index_switches']['Lost'][i+1]=up
        if any('Major' in d['tag'] for d in t):
            indextomajor.append(i+1)
        if any('Minor' in d['tag'] for d in t):
            indextominor.append(i+1)
        if any('Unique' in d['tag'] for d in t):
            indextounique.append(i+1)

    return indexloss,indextomajor,indextominor,indextounique,vis

def analyze_majors(Raw_Results,vis):
    majorloss=[]
    majortoindex=[]
    majortominor=[]
    majortounique=[]

    for i in range (1,len(analyse.results)):
        t=check_raw(Raw_Results,i)
        up=insert_data(t)
        vis["MoSwa_Output"]["Report"]["Results"]['Major_Motifs']['Major_Raw'][i]=up

        if any(not d['hits'] for d in t):
            majorloss.append(i+1)
            z=[d for i,d in enumerate(t) if not d['hits']]
            up=insert_data(z)
            vis["MoSwa_Output"]["Report"]["Results"]['Major_Motifs']['Major_switches']['Lost'][i+1]=up
        if any('Index' in d['tag'] for d in t):
            majortoindex.append(i+1)
        if any('Minor' in d['tag'] for d in t):
            majortominor.append(i+1)
        if any('Unique' in d['tag'] for d in t):
            majortounique.append(i+1)


    return majorloss,majortoindex,majortominor,majortounique,vis

def analyze_minors(Raw_Results,vis):

    minorloss=[]
    minortoindex=[]
    minortomajor=[]
    minortounique=[]

    for i in range (1,len(analyse.results)):
        t=check_raw(Raw_Results,i)
        up=insert_data(t)
        vis["MoSwa_Output"]["Report"]["Results"]['Minor_Motifs']['Minor_Raw'][i]=up

        if any(not d['hits'] for d in t):
            minorloss.append(i+1)
            z=[d for i,d in enumerate(t) if not d['hits']]
            up=insert_data(z)
            vis["MoSwa_Output"]["Report"]["Results"]['Minor_Motifs']['Minor_switches']['Lost'][i+1]=up
        if any('Index' in d['tag'] for d in t):
            minortoindex.append(i+1)
        if any('Major' in d['tag'] for d in t):
            minortomajor.append(i+1)
        if any('Unique' in d['tag'] for d in t):
            minortounique.append(i+1)

    return minorloss,minortoindex,minortomajor,minortounique,vis

def analyze_uniques(Raw_Results,vis):

    uniqueloss=[]
    uniquetoindex=[]
    uniquetomajor=[]
    uniquetominor=[]
    for i in range (1,len(analyse.results)):
        t=check_raw(Raw_Results,i)
        up=insert_data(t)
        vis["MoSwa_Output"]["Report"]["Results"]['Unique_Motifs']['Unique_Raw'][i]=up

        if any(not d['hits'] for d in t):
            uniqueloss.append(i+1)
            z=[d for i,d in enumerate(t) if not d['hits']]
            up=insert_data(z)
            vis["MoSwa_Output"]["Report"]["Results"]['Unique_Motifs']['Unique_switches']['Lost'][i+1]=up
        if any('Index' in d['tag'] for d in t):
            uniquetoindex.append(i+1)
        if any('Major' in d['tag'] for d in t):
            uniquetomajor.append(i+1)
        if any('Minor' in d['tag'] for d in t):
            uniquetominor.append(i+1)
    return uniqueloss,uniquetoindex,uniquetomajor,uniquetominor,vis

def insert_gain(Raw_Results,position,c):
    d={}
    t=check_raw(Raw_Results,position)
    if t:
        t=next((d for i,d in enumerate(t) if c[0] in re.sub('[()]', '', d['seq'])),None)
        d['position']=t['position']
        d['seq']=re.sub('[()]', '', t['seq'])
        d['frequency']=t['count']
        d['tag']=t['tag']
    return d

def gain(motif,vis,Raw_Indexes,Raw_Majors):
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
            if motif == 'Index':
                Raw_Results=Raw_Indexes
                up=insert_gain(Raw_Results,x+1,c)
                vis["MoSwa_Output"]["Report"]["Results"]['Index_Motifs']['Index_switches']['Gain'][x+1]=[up]
            else:
                Raw_Results=Raw_Majors
                up=insert_gain(Raw_Results,x+1,c)
                vis["MoSwa_Output"]["Report"]["Results"]['Major_Motifs']['Major_switches']['Gain'][x+1]=[up]
    return gain,vis

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
                    Raw_Uniques,
                    majortoindex,
                    minortoindex,
                    uniquetoindex,
                    vis):

    mtoi,mitoi,utoi=[],[],[]
    merge_mtoi,merge_mitoi,merge_utoi=[],[],[]

    for i in Filtered_Indexes:

        multi_in = [d for d in Raw_Indexes if d['position'] == i-1 and ' x ' in d['tag']] #get all motif changes related to pos//since we are only filtering only one at a position

        if len(multi_in) < len([d for d in Raw_Indexes if d['position'] == i-1]): #if there are 4 indexes, and 3 of them has switches, no need to look what became an index
            continue

        if i in majortoindex:#if at the same position index is switching to something and
            mtoi.append(i)    #if major is replacing it true switch else it is merging
            t=check_raw(Raw_Majors,i-1)
            z=[d for i,d in enumerate(t) if 'Index' in d['tag']]
            up=insert_data(z)
            vis["MoSwa_Output"]["Report"]["Results"]['Major_Motifs']['Major_switches']['ToIndex'][i]=up
            continue

        if check_raw(Raw_Majors,i-1):
            if any('Index' in str(d['hits'].keys()) for d in check_raw(Raw_Majors,i-1)):
                merge_mtoi.append(i)
                continue

        if i in minortoindex:
            mitoi.append(i)
            t=check_raw(Raw_Minors,i-1)
            z=[d for i,d in enumerate(t) if 'Index' in d['tag']]
            up=insert_data(z)
            vis["MoSwa_Output"]["Report"]["Results"]['Minor_Motifs']['Minor_switches']['ToIndex'][i]=up
            continue

        if check_raw(Raw_Minors,i-1):
            if any('Index' in str(d['hits'].keys()) for d in check_raw(Raw_Minors,i-1)):
                merge_mitoi.append(i)

        if i in uniquetoindex:
            utoi.append(i)
            t=check_raw(Raw_Uniques,i-1)
            z=[d for i,d in enumerate(t) if 'Index' in d['tag']]
            up=insert_data(z)
            vis["MoSwa_Output"]["Report"]["Results"]['Unique_Motifs']['Unique_switches']['ToIndex'][i]=up
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
    merge_mtoi.sort()
    merge_mitoi.sort()
    merge_utoi.sort()


    for x in merge_mtoi:
        t=check_raw(Raw_Majors,x-1)
        z=[d for j,d in enumerate(t) if 'Index' in str(d['hits'].keys())]
        up=insert_data(z)
        vis["MoSwa_Output"]["Report"]["Results"]['Major_Motifs']['Major_Merger']['ToIndex'][x]=up
    for x in merge_mitoi:
        t=check_raw(Raw_Minors,x-1)
        z=[d for j,d in enumerate(t) if 'Index' in str(d['hits'].keys())]
        up=insert_data(z)
        vis["MoSwa_Output"]["Report"]["Results"]['Minor_Motifs']['Minor_Merger']['ToIndex'][x]=up
    for x in merge_utoi:
        t=check_raw(Raw_Uniques,x-1)
        z=[d for j,d in enumerate(t) if 'Index' in str(d['hits'].keys())]
        up=insert_data(z)
        vis["MoSwa_Output"]["Report"]["Results"]['Unique_Motifs']['Unique_Merger']['ToIndex'][x]=up

    return mtoi,mitoi,utoi,merge_mtoi,merge_mitoi,merge_utoi,vis

def m_splits_mergers(Filtered_Majors,
                    Raw_Majors,
                    Raw_Indexes,
                    Raw_Minors,
                    Raw_Uniques,
                    indextomajor,
                    minortomajor,
                    uniquetomajor,
                    vis):

    split_itom,itom,mitom,merge_mitom,utom,merge_utom=([] for i in range(6))
    for i in Filtered_Majors:
        multi_in = [d for d in Raw_Majors if d['position'] == i-1 and ' x ' in d['tag']]

        if len(multi_in) < len([d for d in Raw_Majors if d['position'] == i-1]):
            continue

        if i in indextomajor:
            itom.append(i)
            t=check_raw(Raw_Indexes,i-1)
            z=[d for i,d in enumerate(t) if 'Major' in d['tag']]
            up=insert_data(z)
            vis["MoSwa_Output"]["Report"]["Results"]['Index_Motifs']['Index_switches']['ToMajor'][i]=up
            continue
        if check_raw(Raw_Indexes,i-1):
            if any('Major' in str(d['hits'].keys()) for d in check_raw(Raw_Indexes,i-1)):
                split_itom.append(i)
                t=check_raw(Raw_Indexes,i-1)
                z=[d for j,d in enumerate(t) if 'Major' in str(d['hits'].keys())]
                up=insert_data(z)
                vis["MoSwa_Output"]["Report"]["Results"]['Index_Motifs']['Index_Split']['ToMajor'][i]=up
                continue
        if i in minortomajor:
            mitom.append(i)
            t=check_raw(Raw_Minors,i-1)
            z=[d for i,d in enumerate(t) if 'Major' in d['tag']]
            up=insert_data(z)
            vis["MoSwa_Output"]["Report"]["Results"]['Minor_Motifs']['Minor_switches']['ToMajor'][i]=up
            continue
        if check_raw(Raw_Minors,i-1):
            if any('Major' in str(d['hits'].keys()) for d in check_raw(Raw_Minors,i-1)):
                merge_mitom.append(i)

                continue
        if i in uniquetomajor:
            utom.append(i)
            t=check_raw(Raw_Uniques,i-1)
            z=[d for i,d in enumerate(t) if 'Major' in d['tag']]
            up=insert_data(z)
            vis["MoSwa_Output"]["Report"]["Results"]['Unique_Motifs']['Unique_switches']['ToMajor'][i]=up
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

    merge_mitom.sort()
    merge_utom.sort()

    for x in merge_mitom:
        t=check_raw(Raw_Minors,x-1)
        z=[d for j,d in enumerate(t) if 'Major' in str(d['hits'].keys())]
        up=insert_data(z)
        vis["MoSwa_Output"]["Report"]["Results"]['Minor_Motifs']['Minor_Merger']['ToMajor'][x]=up

    for x in merge_utom:
        t=check_raw(Raw_Uniques,x-1)
        z=[d for j,d in enumerate(t) if 'Major' in str(d['tag'])]
        up=insert_data(z)
        vis["MoSwa_Output"]["Report"]["Results"]['Unique_Motifs']['Unique_Merger']['ToMajor'][x]=up

    return split_itom,itom,mitom,merge_mitom,utom,merge_utom,vis

def mi_splits_mergers(analyse,indextominor,majortominor,uniquetominor,Raw_Indexes,Raw_Majors,Raw_Minors,Raw_Uniques,vis):

    split_itomi,split_mtomi,mtomi,utomi,merge_utomi,migain=([] for i in range(6))
    for x in range (1,len(analyse.results)-1):
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
                        t=check_raw(Raw_Indexes,x)
                        z=[d for i,d in enumerate(t) if 'Minor' in d['tag']]
                        up=insert_data(z)
                        vis["MoSwa_Output"]["Report"]["Results"]['Index_Motifs']['Index_switches']['ToMinor'][x+1]=up
                        continue
                    else:
                        split_itomi.append(x+1)
                        t=check_raw(Raw_Indexes,x)
                        z=[d for j,d in enumerate(t) if 'Minor' in str(d['hits'].keys())]
                        up=insert_data(z)
                        vis["MoSwa_Output"]["Report"]["Results"]['Index_Motifs']['Index_Split']['ToMinor'][x+1]=up
                elif i in a['Major']:
                    if x+1 in majortominor:
                        t=check_raw(Raw_Majors,x)
                        z=[d for i,d in enumerate(t) if 'Minor' in d['tag']]
                        up=insert_data(z)
                        vis["MoSwa_Output"]["Report"]["Results"]['Major_Motifs']['Major_switches']['ToMinor'][x+1]=up
                        continue
                    else:
                        split_mtomi.append(x+1)
                        t=check_raw(Raw_Majors,x)
                        z=[d for j,d in enumerate(t) if 'Minor' in str(d['hits'].keys())]
                        up=insert_data(z)
                        vis["MoSwa_Output"]["Report"]["Results"]['Major_Motifs']['Major_Split']['ToMinor'][x+1]=up
                elif i in a['Unique']:
                    #if x+1 in uniquetominor:
                    utomi.append(x+1)
                    t=check_raw(Raw_Minors,x)
                    z=[d for i,d in enumerate(t) if 'Minor' in d['tag']]
                    up=insert_data(z)
                    vis["MoSwa_Output"]["Report"]["Results"]['Unique_Motifs']['Unique_switches']['ToMinor'][x+1]=up
                else:
                    migain.append(x+1)
                    up=insert_gain(Raw_Minors,x+1,[i])
                    if x+1 not in vis["MoSwa_Output"]["Report"]["Results"]['Minor_Motifs']['Minor_switches']['Gain']:
                        vis["MoSwa_Output"]["Report"]["Results"]['Minor_Motifs']['Minor_switches']['Gain'][x+1]=[]
                    vis["MoSwa_Output"]["Report"]["Results"]['Minor_Motifs']['Minor_switches']['Gain'][x+1].append(up)
    for x in uniquetominor:
        if x not in utomi:
            merge_utomi.append(x)
    utomi=list(set(utomi))
    merge_utomi=list(set(merge_utomi))
    merge_utomi.sort()
    migain=list(set(migain))


    for x in merge_utomi:
        t=check_raw(Raw_Uniques,x-1)
        z=[d for j,d in enumerate(t) if 'Minor' in str(d['tag'])]
        up=insert_data(z)
        vis["MoSwa_Output"]["Report"]["Results"]['Unique_Motifs']['Unique_Merger']['ToMinor'][x]=up


    return split_itomi,split_mtomi,mtomi,utomi,merge_utomi,migain,vis

def u_splits_mergers(analyse,indextounique,majortounique,minortounique,Raw_Indexes,Raw_Majors,Raw_Minors,Raw_Uniques,vis):

    split_itou,split_mtou,split_mitou,ugain=([] for i in range(4))

    for x in range (1,len(analyse.results)-1):
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
                        t=check_raw(Raw_Indexes,x)
                        z=[d for i,d in enumerate(t) if 'Unique' in d['tag']]
                        up=insert_data(z)
                        vis["MoSwa_Output"]["Report"]["Results"]['Index_Motifs']['Index_switches']['ToUnique'][x+1]=up
                        continue
                    else:
                        split_itou.append(x+1)
                        t=check_raw(Raw_Indexes,x)
                        z=[d for j,d in enumerate(t) if 'Unique' in str(d['hits'].keys())]
                        up=insert_data(z)
                        vis["MoSwa_Output"]["Report"]["Results"]['Index_Motifs']['Index_Split']['ToUnique'][x+1]=up
                elif i in a['Major']:
                    if x+1 in majortounique:
                        t=check_raw(Raw_Majors,x)
                        z=[d for i,d in enumerate(t) if 'Unique' in d['tag']]
                        up=insert_data(z)
                        vis["MoSwa_Output"]["Report"]["Results"]['Major_Motifs']['Major_switches']['ToUnique'][x+1]=up
                        continue
                    else:
                        split_mtou.append(x+1)
                        t=check_raw(Raw_Majors,x)
                        z=[d for j,d in enumerate(t) if 'Unique' in str(d['hits'].keys())]
                        up=insert_data(z)
                        vis["MoSwa_Output"]["Report"]["Results"]['Major_Motifs']['Major_Split']['ToUnique'][x+1]=up
                elif i in a['Minor']:
                    if x+1 in minortounique:
                        t=check_raw(Raw_Minors,x)
                        z=[d for i,d in enumerate(t) if 'Unique' in d['tag']]
                        up=insert_data(z)
                        vis["MoSwa_Output"]["Report"]["Results"]['Minor_Motifs']['Minor_switches']['ToUnique'][x+1]=up
                        continue
                    else:
                        split_mitou.append(x+1)
                        t=check_raw(Raw_Minors,x)
                        z=[d for j,d in enumerate(t) if 'Unique' in str(d['hits'].keys())]
                        up=insert_data(z)
                        vis["MoSwa_Output"]["Report"]["Results"]['Minor_Motifs']['Minor_Split']['ToUnique'][x+1]=up
                else:
                    ugain.append(x+1)
                    up=insert_gain(Raw_Uniques,x+1,[i])
                    if x+1 not in vis["MoSwa_Output"]["Report"]["Results"]['Unique_Motifs']['Unique_switches']['Gain']:
                        vis["MoSwa_Output"]["Report"]["Results"]['Unique_Motifs']['Unique_switches']['Gain'][x+1]=[]
                    vis["MoSwa_Output"]["Report"]["Results"]['Unique_Motifs']['Unique_switches']['Gain'][x+1].append(up)

    split_mitou=list(set(split_mitou))
    ugain=list(set(ugain))

    return split_itou,split_mtou,split_mitou,ugain,vis
