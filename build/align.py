import csv,os
from Bio import pairwise2
from Bio.Align.substitution_matrices import load as ld

def align_i(analyse,
        Filtered_Indexes,
        indextomajor,
        indextominor,
        indextounique,
        majortoindex,
        minortoindex,
        uniquetoindex,
        igain,
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
        csvwriter.writerow(['Position', 'Index', 'Motif', 'Score'])

        for x in fix_align:
            tag = ''
            if x in igain:
                tag='Gain'
            elif x in majortoindex:
                tag='Major'
            elif x in minortoindex:
                tag='Minor'
            elif x in uniquetoindex:
                tag='Uniqe'
            csvwriter.writerow([x, fix_align[x]['Index']+' ({})'.format(tag),fix_align[x][list(fix_align[x])[1]]+' ({})'.format(list(fix_align[x])[1]), fix_align[x]['score']])

    return align,fix_align
