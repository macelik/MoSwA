# Motif Switch Analyzer

Takes the multiple sequence alignment of proteins, in fasta, fa, or afa format.


Usage:
```MoSwA.py -i AlignedFasta -m motif_sites -t int -k int -o OutFolderName```

Accepted arguments for -m: index,major,minor,unique (default is all), -t is the support threshold to be analyzed
-k default is 9 for Kmer Length

A summary report is produced at the end.

All the outputs can be found in the created directory.

### Installing Pygraphviz on Windows
1- Download and install 2.46.0 for Windows 10 (64-bit): 
[stable_windows_10_cmake_Release_x64_graphviz-install-2.46.0-win64.exe.](https://gitlab.com/graphviz/graphviz/-/package_files/6164164/download)

2- Install PyGraphviz via
```
python -m pip install --global-option=build_ext --global-option="-IC:\Program Files\Graphviz\include" --global-option="-LC:\Program Files\Graphviz\lib" pygraphviz
```
## Documentation
Protein sequence diversity is one of the major challenges in the design of diagnostic, prophylactic and therapeutic interventions against viruses. Shannon’s entropy has been
used as a quantitative measure of protein sequence diversity, applied via a user-defined k-mer sliding window [1, 2]. The entropy value of a given k-mer position in a sequence
alignment is affected by the total number of distinct k-mer peptides observed at the position and their relative frequency, with a high entropy value indicating diversity,
whereas an entropy of zero indicates complete conservation. Studies [3, 4] have classified distinct k-mer peptides at a given position into diversity motifs (index, major, minor, unique and k-mer_types) based on their incidence (i.e., frequency). Index is the predominant sequence, and all others comprise as total variants, sub-classified into major variant (the predominant variant), minor variants (comprising of k-mers with incidence lower than major and higher than unique) and unique variants (k-mers seen only once in the alignment). These diversity motifs enable diversity dynamics analyses of virus proteins within and between species. Motif switching at a given k-mer alignment position is a phenomenon where fitness change in one or more amino acids, such as through mutations, changes the incidence of a given k-mer sequence across its overlapping positions, resulting in a sequence rank change, and thus, a motif change. This has been observed to occur for all the motifs (index, major, minor or unique). Identifying k-mer positions that exhibited a motif switch and determining the nature of the switches was a challenge given the large combination of switches that are possible and their omnipresence. Herein, we present MoSwA (Motif Switch Analyzer; https://github.com/macelik/MoSwA), a tool that not only identifies all alignment k-mer positions that exhibit motif switching, but also provides a multi-faceted and extensive characterisation of the switches. This includes:
-	a short statistical summary on the switches observed;
-	an alignment view of all the switches observed in a given dataset, referenced against a consensus sequence built from the index sequences of the k-mer positions;
-	a network graph showing the dynamic interaction between the motif switches;
-	a bar plot showing clusters of motif switch positions, as well as hotspots in the protein alignment;
-	index switch positions are noteworthy and highlighted because highly conserved index at such positions are possibly to be avoided as vaccine targets given the instability of the index;
-	a pairwise alignment score based on PAM30 is provided for index switches, to determine the physico-chemical spectrum of similarity/variability between the index sequence and the replacing variant motif sequence. 
The input to MoSwA is a protein multiple sequence alignment. The user can choose to study a single, multiple or all diversity motifs. MoSwA now enables a comparative analyses of motif switches within and between viral species proteomes. This can help provide a better understanding of motif switches in relation to viral sequence diversity and evolution, with implications to vaccine and drug design.


### Sample Output

![Cluster_and_Hotspots](/showface/Allin1.PNG)
Figure 1: Selected features of MoSwA. A) A short statistical summary on the switches observed. B) A bar plot showing clusters of motif switch positions, as well as hotspots in the protein alignment. C) A network graph showing the dynamic interaction between the motif switches.


### Sample of dynamic network dot plot

[Please go to the link to view the example dynamic plot](/showface/Dynamic_Plot.html)





