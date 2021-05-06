# Motif Switch Analyzer

Takes the json output from Hunana.


Usage:
```motifswitcher.py -i jsonFileInput -m motif_sites -t int```

Accepted arguments for -m: index,major,minor,unique (default is all), -t is the support threshold to be analyzed

A summary report is produced at the end.
Image.png and Output.fasta must be in the same directory with the report.


## Documentation

  Takes the input of AVANA or HUNANA, first strip down the json file to simplify it in order to hold less storage.
  A given position has one Index and Major but might have more than one Minor and Unique. Thus number tags are needed for these motifs. We tag them with an incrementing number, the one with the highest frequency will be tagged as Minor_1 and randomly for Unique since Unique motifs have frequency of only 1.

### Find overlapping motifs
  Then we seperate the specific motifs at a given position which allows us to compare seqs, tags, counts and positions. We drop the first amino acid in the preceding position while dropping the last amino acid in the next position, Each motif in the preceding position treated seperately and compared to everything in the next position. The output is procudes for each motif in a JSON format to the user.

### Identify motif changes

  If a motif tag of index and major of a current position is not found in the next position with the same seq, we mark it as a motif switch. However for minor motifs, we not only check the sequence and the tag but also frequency.  To tag a minor as a motif switch there has to be tag and count difference since minor tag can change but frequency and seq might stay the same. or it has to come out out of minor tag then regardless of the count, we identify it as a switch.

  Later, these indentified motif switch positions are filtered and given to the user as a JSON file which includes preciding position and next position, seqs, frequencies and motifs that are switched.


### Concataneting A Consensus Sequence from The Index:

  A reverse sliding (k-mer size of the motif-1) window approach was utilized to identify overlapping amino acid at a position. To identify a consensus amino acid at a position, indicies of preciding k-mer positions are taken with their frequencies. The amino acid with the highest frequency is determined as consensus a.a at that position. If overlapping a.a for that position is %100, shown as capital letter otherwise shown as small letter.
```
          ^
          |
  GLTGTTTVT			pos=178
   LTGTTTVTP
    NGTNGTNQE
     GTNGTNQER
      TNGTNQERT
       NGTNQERTN
        GTNQERTNT
         TNQERTNTT
          NQERTNTTN	pos=170

```

Along with built index consensus, we also show the motif changes in a Clustalw format as below.
Asterix(*) indicates motif switch, while colon(:) indicates motif switch position with a low support value.



```
index_consensus                     MAYSSSRLLIALLLiSiygivcatkpkrqyvtvfYGIPAWRNASIPLFCA
index_switches                      MAYSSSRLLIALL*iSiy:ivcatkpkrqyvtvfYGIPAWRNASIPLFCA
major_switches                      ***S**R*L**L**iSiygivcatkpkr****vfYGIP*W*NASIPLFCA
minor_switches                      **************iSiygivcatkpkr****vf***P*W*N*S*PLFCA
unique_switches                     MA****RLL**L**iSiygivcatkpkr*yv*vfYGIP*W*NA*IPL*CA

index_consensus                     TKNRDTWGTIQCLPDNDDYQEIALNVTEAFDAWdNTVTEQAVEDVWSLFE
index_switches                      TKNRDTWGTIQCLPDNDDYQEIALN**EAFDAWdNTVTEQAVEDVWSLFE
major_switches                      TKNRDTWGTIQCL*DND**QEI*L***EAFDAWdNTVTEQAVEDVW*LFE
minor_switches                      T*NRDTW*T**CL**N****E***N*TEAFDA**N**TEQA**DVW*LFE
unique_switches                     T*NRDTWG***CLPDN**Y*EI***V*EAFDAW**TVTEQA*EDV*SLFE
```

### Clusters and Hotspots

  Unique motif switch positions are taken to indetify clusters. We check the distances between kmers to form a cluster, distance is determined by subtraction of preciding and next position. If distance is equal or less than k-mer motif size, the positions are placed in the same cluster. Another criteria to form a cluster is that a cluster must have more than the minimum number of switch kmers (5).

  To form hotspot, the clusters must be within a certain distance and there has to be at least 2 clusters.
  Cmsp; is the number of positions that make up a cluster while Cs is the cluster size. Cd is the distance between clusters, The starting and ending positions that make up the cluster shown in square brackets as [start::end]. Red square brackets point out to motif switch hotspots.

  ![Cluster_and_Hotspots](/showface/clusters_hotpost.png)


### Plotting Motif Switches

After identifying motif swithces, the question becomes, what motif replaced the it? We indetify these too and plot them dynamically.
[Please go to the link to view the example dynamic plot](/showface/Dynamic_Plot.html)





