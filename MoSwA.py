#!/usr/bin/env python
from argparse import ArgumentParser
from build import RunAndCompare,SplitsMergers,BuildConsensus,Clustering,PlotHotClusters
from build import NetworkPlot,report,align
import os,random,json

def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg  # return an open file handle

def check_positive(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise parser.error("%s is an invalid positive int value" % value)
    return ivalue

def create_dir(out):
    path = os.path.join(os.getcwd(), out)
    try:
        os.mkdir(path)
        print("Directory {} created successfully".format(out))
    except OSError as error:
        print("Directory '%s' exists" % out)
        rand=random.randint(1000,5000)
        path=path+"_"+str(rand)
        os.mkdir(path)
        print("Directory {}_{} created instead" .format(out,rand))

    return path



parser = ArgumentParser(description="Motif finder")
parser.add_argument("-i", dest="filename", required=True,
                    help="json file input", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))
#parser.add_argument("-k", default=9, type=int, help="This is the k-mer size")
parser.add_argument("-m",
                    choices=["index", "major", "minor", "unique","all"],
                    nargs="+",
                    default="all", type=str, help="motif sites to analyze")
parser.add_argument("-t", type=check_positive,default=30,
                    help="Minimum number of seqeunces at a position")
parser.add_argument("-k", type=check_positive,default=9,
                    help="Minimum number of seqeunces at a position")
parser.add_argument("-o", help="Output directory name", required=True)

def create_dic():
    dict1={}
    dict1['I => M']=indextomajor+split_itom
    dict1['M => I']=mtoi+merge_mtoi
    dict1['I => Mi']=indextominor+split_itomi
    dict1['Mi => I']=mitoi+merge_mitoi
    dict1['I => U']=indextounique+split_itou
    dict1['U => I']=utoi+merge_utoi
    dict1['I => +']=igain
    dict1['I => -']=indexloss
    dict1['M => Mi']=majortominor+split_mtomi
    dict1['Mi => M']=mitom+merge_mitom
    dict1['M => U']=majortounique+split_mtou
    dict1['U => M']=utom+merge_utom
    dict1['M => +']=mgain
    dict1['M => -']=majorloss
    dict1['Mi => U']=minortounique+split_mitou
    dict1['U => Mi']=utomi+merge_utomi
    dict1['Mi => +']=migain
    dict1['Mi => -']=minorloss
    dict1['U => +']=ugain
    dict1['U => -']=uniqueloss
    return dict1

args = parser.parse_args()
dima_kmer = args.k
which=args.m
thold=args.t
fname=args.filename
out=args.o
path=create_dir(out)

if 'all' in which:
    which='index','major','minor','unique'


RunAndCompare.path = path
analyse = RunAndCompare.GetMotifs(fname,dima_kmer)
common_pos,lengthofswitches,Raw_Indexes,Filtered_Indexes,Raw_Majors,Filtered_Majors,Raw_Minors,Filtered_Minors,Raw_Uniques,Filtered_Uniques=RunAndCompare.runit(which,dima_kmer,analyse)
common_pos.sort()
SplitsMergers.analyse = analyse
vis=SplitsMergers.view_json()
indexloss,indextomajor,indextominor,indextounique,vis = SplitsMergers.analyze_index(Raw_Indexes,vis)
majorloss,majortoindex,majortominor,majortounique,vis = SplitsMergers.analyze_majors(Raw_Majors,vis)
minorloss,minortoindex,minortomajor,minortounique,vis = SplitsMergers.analyze_minors(Raw_Minors,vis)
uniqueloss,uniquetoindex,uniquetomajor,uniquetominor,vis = SplitsMergers.analyze_uniques(Raw_Uniques,vis)
mtoi,mitoi,utoi,merge_mtoi,merge_mitoi,merge_utoi,vis = SplitsMergers.i_splits_mergers(Filtered_Indexes,Raw_Indexes,Raw_Majors,Raw_Minors,Raw_Uniques,majortoindex,minortoindex,uniquetoindex,vis)
igain,vis = SplitsMergers.gain('Index',vis,Raw_Indexes,Raw_Majors)
mgain,vis = SplitsMergers.gain('Major',vis,Raw_Indexes,Raw_Majors)
split_itom,itom,mitom,merge_mitom,utom,merge_utom,vis = SplitsMergers.m_splits_mergers(Filtered_Majors,
                                                                                    Raw_Majors,
                                                                                    Raw_Indexes,
                                                                                    Raw_Minors,
                                                                                    Raw_Uniques,
                                                                                    indextomajor,
                                                                                    minortomajor,
                                                                                    uniquetomajor,
                                                                                    vis)
split_itomi,split_mtomi,mtomi,utomi,merge_utomi,migain,vis = SplitsMergers.mi_splits_mergers(analyse,indextominor,majortominor,uniquetominor,Raw_Indexes,Raw_Majors,Raw_Minors,Raw_Uniques,vis)
split_itou,split_mtou,split_mitou,ugain,vis = SplitsMergers.u_splits_mergers(analyse,indextounique,majortounique,minortounique,Raw_Indexes,Raw_Majors,Raw_Minors,Raw_Uniques,vis)
to_build,no_sup,low_sup,highest,average = BuildConsensus.building_index(analyse,thold,dima_kmer)
cons = BuildConsensus.ConsensusInput.parse_consensus_input(to_build)
cons = BuildConsensus.calculate_consensus(consensus_input=cons)
unable_to_analyze = BuildConsensus.unable_to_analyze(analyse,Filtered_Indexes,Filtered_Majors,dima_kmer,cons)
Filtered_Uniques,Filtered_Minors,Filtered_Majors,Filtered_Indexes = BuildConsensus.check_which(which,Filtered_Indexes,Filtered_Majors,Filtered_Minors,Filtered_Uniques)
consensus = BuildConsensus.build_consensus(Filtered_Indexes, Filtered_Majors, Filtered_Minors, Filtered_Uniques, to_build, thold, os.path.join(path, "Consensus.aln"))
cluster,r_cluster = Clustering.cluster_it(common_pos,no_sup,low_sup)
Clustering.CnandCs(r_cluster,dima_kmer)
print("Plotting Clusters and Hotspots")
hots,r_hots = Clustering.hotspots(r_cluster)
plotin = PlotHotClusters.plotit(r_cluster,dima_kmer)
fig = PlotHotClusters.draw(r_cluster,r_hots,plotin,path)
print("Plotting Network Graph")
dict1 = create_dic()
dict1 = NetworkPlot.check_arg(dict1,which)
repos,D,connection = NetworkPlot.network(dict1)

sandm={'split_itom':split_itom,'split_itomi':split_itomi,'split_itou':split_itou,
       'split_mtomi':split_mtomi,'split_mtou':split_mtou,'split_mitou':split_mitou,
       'merge_mtoi':merge_mtoi,'merge_mitoi':merge_mitoi,'merge_utoi':merge_utoi,
       'merge_mitom':merge_mitom,'merge_utom':merge_utom,'merge_utomi':merge_utomi}

s_edges,s_nodes = NetworkPlot.plotly_edges(repos,D,connection,sandm)
data = NetworkPlot.plotly_info(s_edges,s_nodes)
NetworkPlot.out_plotly(data,path)
NetworkPlot.fix_html(path)
alignment_length=len(analyse.results)+dima_kmer-1
summed=report.write_report(average,thold,no_sup,low_sup,common_pos,lengthofswitches,unable_to_analyze,alignment_length,highest,path)
vis=report.report_json(average,thold,no_sup,low_sup,common_pos,lengthofswitches,unable_to_analyze,alignment_length,highest,summed,dima_kmer,vis)
align,fix_align,vis = align.align_i(analyse,Filtered_Indexes,indextomajor,merge_mtoi,indextominor,merge_mitoi,indextounique,majortoindex,minortoindex,uniquetoindex,igain,vis,path)

if "index" not in which:
    del vis['MoSwa_Output']['Report']['Results']['Index_Motifs']
if "major" not in which:
    del vis['MoSwa_Output']['Report']['Results']['Major_Motifs']
if "minor" not in which:
    del vis['MoSwa_Output']['Report']['Results']['Minor_Motifs']
if "unique" not in which:
    del vis['MoSwa_Output']['Report']['Results']['Unique_Motifs']

f = open(os.path.join(path, "Results.json"),'w')
json.dump(vis,f)
f.close()
