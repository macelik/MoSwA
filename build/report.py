import sys
def write_report(myd,average,no_sup,low_sup,thold,topla,indexes,majors,minors,uniques,indextomajor,majortoindex,indextominor,minortoindex,
                indextounique,uniquetoindex,indexgain,indexloss,majortominor,minortomajor,majortounique,uniquetomajor,majorgain,majorloss,minortounique,
                uniquetominor,minorgain,minorloss,uniquegain,uniqueloss,common_pos,unable_to_analyze):
    highest_sup=0
    for x in range(0,len(myd)):
        if myd[x]['supports'] > highest_sup:
            highest_sup=myd[x]['supports']

    if indextounique==[]:
        indextounique='-'
    if unable_to_analyze ==[]:
        unable_to_analyze='NA'




    replacements = {'alignment_length!':str(myd[-1]['position']), 'longest_support!':str(highest_sup),
                'average!':str(average), 'len_no_sup!':str(len(no_sup)),'len_low_sup!':str(len(low_sup)),
                '!thold!':str(thold),'uniqueswitches!':str(len(topla)),'len_indexes!':str(len(indexes)),
                'len_majors!':str(len(majors)),'len_minors!':str(len(minors)),
                'len_uniques!':str(len(uniques)),'IndextoMajor':', '.join(map(str,indextomajor)),
                'MajortoIndex':', '.join(map(str,majortoindex)),
                'IndextoMinor':', '.join(map(str,indextominor)),'MinortoIndex':', '.join(map(str,minortoindex)),
                'IndextoUnique':', '.join(map(str,indextounique)),'UniquetoIndex':', '.join(map(str,uniquetoindex)),
                'IndexGain':', '.join(map(str,indexgain)),'IndexLoss':', '.join(map(str,indexloss)),
                'MajortoMinor':', '.join(map(str,majortominor)),'MinortoMajor':', '.join(map(str,minortomajor)),
                'MajortoUnique':', '.join(map(str,majortounique)),'UniquetoMajor':', '.join(map(str,uniquetomajor)),'MajorGain':', '.join(map(str,majorgain)),
                'MajorLoss':', '.join(map(str,majorloss)),'MinortoUnique':', '.join(map(str,minortounique)),
                'UniquetoMinor':', '.join(map(str,uniquetominor)),'MinorGain':', '.join(map(str,minorgain)),'MinorLoss':', '.join(map(str,minorloss)),
                'UniqueGain':', '.join(map(str,uniquegain)),'UniqueLoss':', '.join(map(str,uniqueloss)),'totalswitches!':str(len(common_pos)),
                 '!UnableToAnalyze':'/'.join(map(str,unable_to_analyze))}

    with open('build/htmltemplate/index.html') as infile, open('report.html', 'w') as outfile:
        for line in infile:
            for src, target in replacements.items():
                line = line.replace(src, target)
            outfile.write(line)