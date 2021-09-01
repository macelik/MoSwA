import json,os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import itertools

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

def plotit(r_cluster,kmer):
    plotin={}
    a=list(r_cluster.keys())
    for r in r_cluster:
        get_index_pos=a.index(r)
        if get_index_pos + 1 == len(a):
            Cs=r_cluster[r][-1][-1] - r_cluster[r][0][0] + kmer
            Cmsp=len(r_cluster[r])+1
            start,end=r_cluster[r][0][0],r_cluster[r][-1][-1]
            plotin[r]={'Cmsp':Cmsp,'Cs':Cs,'Dist':None,'Start':start,'End':end}

        else:
            get_next_cluster=a[get_index_pos+1]
            Cs=r_cluster[r][-1][-1] - r_cluster[r][0][0] + kmer
            Cmsp=len(r_cluster[r])+1
            dist= abs(r_cluster[get_next_cluster][0][0] - r_cluster[r][-1][-1])
            start,end=r_cluster[r][0][0],r_cluster[r][-1][-1]
            plotin[r]={'Cmsp':Cmsp,'Cs':Cs,'Dist':dist,'Start':start,'End':end}
    return plotin

def draw(r_cluster,r_hots,plotin,path):
    get_pos_cluster={}
    fig, ax = plt.subplots(figsize=(18, 11))
    hc_lst =([plotin[x]['Cs'] for x in plotin])
    a=[plotin[x]['Cmsp'] for x in plotin]
    index=list(range(0,len(plotin)*3,3))
    bar_width = 1
    rects1 = plt.bar(index, hc_lst, len(r_cluster)/10, color='g',label='Cluster Size (Cs)', alpha= 0.5, align='center')
    rects2 = plt.bar(index, a, len(r_cluster)/10, color='g',label='Number of Positions (Cmsp)', alpha= 0.8, align='center')
    rects3 = plt.bar(index, 0, color='w', align='center')

    plt.ylabel('Cluster Size in Bp')
    plt.xlabel('Number of Clusters')

    plt.xticks(index, ([key for key in plotin]))
    plt.legend()
    for myBar,keys in zip(rects1,plotin):
        get_pos_cluster[keys]=myBar.get_x()
        midleTextX = myBar.get_x() + bar_width
        midleTextY = 2
        ax.text(myBar.get_x() + myBar.get_width()/2.0 ,myBar.get_height(), 'Cmsp: {} \n[{}::{}]\nCs: {}'.format(plotin[keys]['Cmsp'],plotin[keys]['Start'],plotin[keys]['End'],plotin[keys]['Cs']),ha='center', va='bottom', rotation='horizontal', ma='left')
        if plotin[keys]['Dist'] == None:
            continue
        else:
            ax.text(midleTextX+(bar_width/2),midleTextY,'Cd: {}\n<==>'.format(plotin[keys]['Dist']),ha='left', va='bottom', rotation='horizontal', ma='left')

    for keys in r_hots:
        start=get_pos_cluster[r_hots[keys][0][0]]-(myBar.get_width()/5)
        end=get_pos_cluster[r_hots[keys][-1][-1]]+myBar.get_width()+(myBar.get_width()/5)
        lst=list(set(itertools.chain(*r_hots[keys])))
        hei=0
        for x in lst:
            if hei < plotin[x]['Cs']:
                hei=plotin[x]['Cs']

        ax.add_patch(
        patches.Rectangle(
            xy=(start,0),  # point of origin.
            width=end-start,
            height=hei+5,
            linewidth=2,
            color='red',
            fill=False
        )
    )
    fig.tight_layout()
    fig.savefig(os.path.join(path, "ClusterHotSpotPlot.png"),bbox_inches='tight',dpi=300)
    return fig
