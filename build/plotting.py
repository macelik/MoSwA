import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import itertools
import datetime
from matplotlib.backends.backend_pdf import PdfPages

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

def draw(r_hots,plotin):
    get_pos_cluster={}
    fig, ax = plt.subplots(figsize=(18, 10))
    hc_lst =([plotin[x]['Cs'] for x in plotin])
    a=[plotin[x]['Cmsp'] for x in plotin]
    index=list(range(0,len(plotin)*3,3))
    bar_width = 1
    rects1 = plt.bar(index, hc_lst, bar_width, color='g',label='Cluster Size (Cs)', alpha= 0.5, align='center')
    rects2 = plt.bar(index, a, bar_width, color='g',label='Number of Positions (Cmsp)', alpha= 0.8, align='center')

    plt.ylabel('Cluster Size in Bp')
    plt.xlabel('Cluster Number')

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
        start=get_pos_cluster[r_hots[keys][0][0]]
        end=get_pos_cluster[r_hots[keys][-1][-1]]
        lst=list(set(itertools.chain(*r_hots[keys])))
        hei=0
        for x in lst:
            if hei < plotin[x]['Cs']:
                hei=plotin[x]['Cs']

        ax.add_patch(
        patches.Rectangle(
            xy=(start,0),  # point of origin.
            width=end-start+1,
            height=hei+12,
            linewidth=2,
            color='red',
            fill=False
        )
    )
    fig.tight_layout()
    fig.savefig("image2.png",bbox_inches='tight',dpi=100)
    return fig

def firstpage(text1,text2,text3):
    fig1, ax1 = plt.subplots(figsize=(11.69,8.27))
    ax1.set_axis_off()
    ax1.set_title('Motif Switch Analyzer Report',fontsize=30, fontweight='bold')

    #text1 =  ""
    #text2 = "Text2"


    ax1.text(-0.12, 0.82,text1, transform=ax1.transAxes, fontsize=12)
    ax1.text(-0.12, 0.68,text2, transform=ax1.transAxes, fontweight='bold')
    ax1.text(-0.12, 0.18,text3, transform=ax1.transAxes)
    return fig1

def write_pdf(fig,fig1):    
    pp = PdfPages('multipagFinal.pdf')
    fig1.savefig(pp, format='pdf', orientation='portrait')
    fig.savefig(pp, format='pdf', orientation='portrait')
    d = pp.infodict()
    d['Title'] = 'Motif Switch Report'
    d['Author'] = 'Mohammad Asif Khan & Muhammet Celik'
    d['Subject'] = 'Finding Motif Switches from Hunana Output'
    d['Keywords'] = 'Some keywords here'
    d['CreationDate'] = datetime.datetime(2009, 11, 13)
    d['ModDate'] = datetime.datetime.today()
    pp.close()