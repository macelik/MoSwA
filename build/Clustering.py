def cluster_it(common_pos,no_sup,low_sup):

    topla=[]

    for x in common_pos:
        if x in no_sup+low_sup:
            continue
        else:
            topla.append(x)
    topla.sort()


    cluster={}
    r_cluster={}
    c_num=1
    cluster['c_'+str(c_num)]=[]
    for i in range (0,len(topla)-1):
        pos1=topla[i]
        pos2=topla[i+1]
        Cd=abs(pos2 - pos1)
        if Cd <= 9:
            cluster['c_'+str(c_num)].append((pos1,pos2))
        if Cd > 9:
            c_num = c_num+1
            cluster['c_'+str(c_num)]=[]
    c_num=0
    for keys in cluster:
        if len(cluster[keys]) >= 6: #this 6 actually means 7 cause counting in touples. remember to add +1
            c_num=c_num+1
            #print(c_num)
            r_cluster['c_'+str(c_num)]=cluster[keys]
            #r_cluster['c_'+str(c_num)].append(cluster[keys])
        else:
            continue
    return cluster,r_cluster

def CnandCs(r_cluster,kmer):
    print('The number of clusters found: {}\nCmsp; is the number of positions that make up a cluster while Cs is the cluster size'.format(len(r_cluster)))
    for r in r_cluster:
        Cs=r_cluster[r][-1][-1] - r_cluster[r][0][0] + kmer
        print ('Cmsp of cluster {} is {} and Cs is'.format(r,len(r_cluster[r])+1),Cs)

def hotspots(r_cluster):
    a=list(r_cluster.keys())
    hots={}
    c_num=1
    hots['hotspot_'+str(c_num)]=[]
    for i in range(0,len(a)-1):
        pos1=r_cluster[a[i]][-1][-1]
        pos2=r_cluster[a[i+1]][0][0]
        Hd=pos2 - pos1
        if Hd < 20:
            hots['hotspot_'+str(c_num)].append((a[i],a[i+1]))
        if Hd > 20:
            c_num=c_num+1
            hots['hotspot_'+str(c_num)]=[]
    c_num=0
    r_hots={}
    for keys in hots:
        if len(hots[keys]) > 0:
            c_num=c_num+1
            #print(c_num)
            r_hots['hotspot_'+str(c_num)]=hots[keys]
            #r_cluster['c_'+str(c_num)].append(cluster[keys])
        else:
            continue
    return hots,r_hots
