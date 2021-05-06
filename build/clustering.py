def low_no_sup(myd,thold):
	low_no_sup=[]
	low_sup=[]
	no_sup=[]
	for x in range(0,len(myd)):
	    if myd[x]['supports'] == 0:
	        no_sup.append(myd[x]['position'])
	        low_no_sup.append(myd[x]['position'])
	        continue
	    if myd[x]['supports'] < thold:
	        low_sup.append(myd[x]['position'])
	        low_no_sup.append(myd[x]['position'])
	return low_no_sup,no_sup,low_sup

def uniqpos_no_low_support(common_pos,low_no_sup): #unique positions low or no sup positions are extracted
	topla=[]
	common_pos=list(set(common_pos))
	for i in common_pos:
	    if i in low_no_sup:
	        continue
	    else:
	        topla.append(i)
	topla.sort()
	return topla

def cluster_num():      
	global c_num
	pStart = 1 #adjust start value, if req'd   
	pInterval = 1 #adjust interval value, if req'd  
	if (c_num == 0):   
		c_num = pStart   
	else:   
		c_num = c_num + pInterval   
	return c_num

def cluster_it(topla):
    global c_num
    cluster={}
    r_cluster={}
    c_num=0
    c_num=cluster_num()
    cluster['c_'+str(c_num)]=[]
    for i in range (0,len(topla)-1):
        pos1=topla[i]
        pos2=topla[i+1]
        Cd=abs(pos2 - pos1)
        if Cd <= 9:
            cluster['c_'+str(c_num)].append((pos1,pos2))
        if Cd > 9:
            c_num=cluster_num()
            cluster['c_'+str(c_num)]=[]
    c_num=0
    for keys in cluster:
        if len(cluster[keys]) >= 6: #this 6 actually means 7 cause counting in touples. remember to add +1
            c_num=cluster_num()
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
    global c_num
    a=list(r_cluster.keys())
    hots={}
    c_num=0
    c_num=cluster_num()
    hots['hotspot_'+str(c_num)]=[]
    for i in range(0,len(a)-1):
        pos1=r_cluster[a[i]][-1][-1]
        pos2=r_cluster[a[i+1]][0][0]
        Hd=pos2 - pos1
        if Hd < 20:
            hots['hotspot_'+str(c_num)].append((a[i],a[i+1]))
        if Hd > 20:
            c_num=cluster_num()
            hots['hotspot_'+str(c_num)]=[]
    c_num=0
    r_hots={}
    for keys in hots:
        if len(hots[keys]) > 0: #this 6 actually means 7 cause counting in touples. remember to add +1
            c_num=cluster_num()
            #print(c_num)
            r_hots['hotspot_'+str(c_num)]=hots[keys]
            #r_cluster['c_'+str(c_num)].append(cluster[keys])
        else:
            continue
    return hots,r_hots