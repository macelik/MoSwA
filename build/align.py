import csv,os
from Bio import pairwise2
from Bio.Align.substitution_matrices import load as ld

def align_i(analyse,
        Filtered_Indexes,
        indextomajor,
        merge_mtoi,
        indextominor,
        merge_mitoi,
        indextounique,
        majortoindex,
        minortoindex,
        uniquetoindex,
        igain,
        vis,
        path):
    matrix = ld(("PAM30"))
    gap_open = -10
    gap_extend = -0.5
    align={}

    for x in indextomajor:
        align[x]={'Index':list(analyse.get_index('Index',x-1).keys())[0][0],'Major':list(analyse.get_index('Major',x-1).keys())[0][0]}
        a=align[x]['Index']
        b=align[x]['Major']
        alns = pairwise2.align.globalds(a, b, matrix, gap_open, gap_extend)
        align[x]['score']=alns[0][2]

    for x in indextominor:
        pos=x
        tag=list(Filtered_Indexes[x].values())[3].split(" ")[-1]
        seq1=list(analyse.get_index('Index',x-1).keys())[0][0]
        for i in list(analyse.get_index('Minor',x-1).keys()):
            if i[1] == tag:
                seq2 = i[0]
                align[x]={'Index':seq1,str(tag):seq2}
                a=align[x]['Index']
                b=align[x][str(tag)]

                alns = pairwise2.align.globalds(a, b, matrix, gap_open, gap_extend)
                align[x]['score']=alns[0][2]

    for x in indextounique:
        pos=x
        tag=list(Filtered_Indexes[x].values())[3].split(" ")[-1]
        if analyse.get_index('Index',x-1):
            seq1=list(analyse.get_index('Index',x-1).keys())[0][0]
            for i in list(analyse.get_index('Unique',x-1).keys()):
                if i[1] == tag:
                    seq2 = i[0]
                    align[x]={'Index':seq1,str(tag):seq2}
                    a=align[x]['Index']
                    b=align[x][str(tag)]

                    alns = pairwise2.align.globalds(a, b, matrix, gap_open, gap_extend)
                    align[x]['score']=alns[0][2]
    align=dict(sorted(align.items()))

    fix_align={}

    for pos in align:
        a=''
        b=''
        for x,y in zip(list(align[pos].values())[0],list(align[pos].values())[1]):
            if x==y:
                a =a+'.'
                b =b+'.'
            else:
                a +=x
                b +=y
        fix_align[pos]={'Index':a,list(align[pos].keys())[1]:b,'score':align[pos]['score'] }

    with open(os.path.join(path, "IndexSwitchAlignment.csv"), "w") as csv_file:
        csvwriter = csv.writer(csv_file)
        csvwriter.writerow(['Position', 'Current Index', 'Motif', 'Score','Support','Index %','Major %','Minor %','Unique %'])

        for x in fix_align:
            tag = ''
            pos_sup=analyse.results[x-1].support
            index_p=list(analyse.get_index('Index',x-1).values())[0][1]*len(analyse.get_index('Index',x-1))
            if len(analyse.get_index('Major',x-1)) == 0:
                major_p=0
            else:
                major_p=list(analyse.get_index('Major',x-1).values())[0][1]*len(analyse.get_index('Major',x-1))
            uniq_p=len(analyse.get_index('Unique',x-1))
            minor_p=pos_sup-(index_p+major_p+uniq_p)
            index_p=round(index_p/pos_sup*100)
            major_p=round(major_p/pos_sup*100)
            minor_p=round(minor_p/pos_sup*100)
            uniq_p=round(uniq_p/pos_sup*100)
            if x in igain:
                tag='Gain'
            elif x in majortoindex or x in merge_mtoi:
                tag='Major'
            elif x in minortoindex or x in merge_mitoi:
                tag='Minor'
            elif x in uniquetoindex:
                tag='Uniqe'
            csvwriter.writerow([x, fix_align[x]['Index']+' ({})'.format(tag),fix_align[x][list(fix_align[x])[1]]+' ({})'.format(list(fix_align[x])[1]), fix_align[x]['score'],pos_sup,index_p,major_p,minor_p,uniq_p])
            a={
            "Index_replaced_with":tag,
            "Index_switched_to":list(fix_align[x])[1],
            "Current_Index":fix_align[x]['Index'],
            "Current_"+str(list(fix_align[x])[1]):fix_align[x][list(fix_align[x])[1]],
            "Score_PAM30":fix_align[x]['score'],
            "Support":pos_sup,
            "Index_%":index_p,
            "Major_%":major_p,
            "Minor_%":minor_p,
            "Unique_%":uniq_p}
            vis['MoSwa_Output']['Report']['PairWise_Alignments'][x]="Null"
            vis['MoSwa_Output']['Report']['PairWise_Alignments'][x]=a
            #print(fix_align)

    return align,fix_align,vis
