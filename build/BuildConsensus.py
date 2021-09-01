from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
from collections import defaultdict


class SwitchData:
    def __init__(self, id_, positions):
        self.id_ = id_
        self.positions = positions

    @staticmethod
    def parse_filtered_data(filtered_data_type: str, data: dict):
        switches = []
        if data == {}:
            switches.append(int(0))
        for key in data:
            switches.append(int(key))

        switch_data = SwitchData("{}_switches".format(filtered_data_type), switches)

        return switch_data


class ConsensusInput:
    def __init__(self, position: int, data: str, support_count: int, count: int):
        self.position = position
        self.data = data
        self.support_count = support_count
        self.count = count

    @staticmethod
    def parse_consensus_input(data):
        output = []
        for value in data:
            output.append(ConsensusInput(value[0], value[1], value[2], value[3]))

        return output


def calculate_consensus(consensus_input: list):
    output = StringIO()
    length_of_each_line = len(consensus_input[0].data)
    total_unique_positions = consensus_input[len(consensus_input) - 1].position - consensus_input[0].position + 1
    for i in range(length_of_each_line + total_unique_positions - 1):
        count = defaultdict(int)
        for j in range(max(0, i - length_of_each_line + 1), len(consensus_input)):
            target_position = consensus_input[j].position - consensus_input[0].position
            if i - target_position >= length_of_each_line:
                continue
            if target_position > i:
                break
            #print(i, target_position, j)
            alphabet = consensus_input[j].data[i - target_position]
            count[alphabet] += consensus_input[j].count
        max_count = 0
        max_alphabet = "-"
        for alphabet in count:
            # Ignore counting for hyphen
            if alphabet == '-':
                continue
            new_count = count[alphabet]
            if new_count > max_count:
                max_alphabet = alphabet
                max_count = new_count
        if len(count) == 1:
            output.write(max_alphabet)
        else:
            output.write(max_alphabet.lower())

    return output.getvalue()


def get_data_from_consensus_file(filepath) -> list:
    output = []
    with open(filepath) as f:
        for line in f:
            output.append(line.strip())
    return output


def write_consensus(filepath,
                 consensus_input,
                 threshold,
                 consensus_output,
                 index_switches,
                 major_switches,
                 minor_switches,
                 unique_switches):
    with open(filepath, 'w') as f:
        records = []
        consensus_record = SeqRecord(
            seq=Seq(consensus_output),
            id="index_consensus",
            description="",
        )
        records.append(consensus_record)

        for switches_data in [index_switches, major_switches, minor_switches, unique_switches]:
            new_consensus_output = list(consensus_output)

            if 0 in switches_data.positions:
                new_consensus_output = 'NA'
                continue

            else:

                for position in switches_data.positions:


                    if position - 1 > len(consensus_output) or position - 1 < 0:
                        continue
                    alphabet_output = "*"
                    if position - 1 < len(consensus_input):
                        if consensus_input[position - 1].support_count < threshold:
                            alphabet_output = ":"
                    new_consensus_output[position - 1] = alphabet_output
                new_consensus_output = ''.join(new_consensus_output)





            record = SeqRecord(
                seq=Seq(new_consensus_output),
                id=switches_data.id_,
                description="",
            )
            records.append(record)

        SeqIO.write(records, filepath, "clustal")


def build_consensus(index_switches, major_switches, minor_switches, unique_switches, to_build, threshold, output_file_name):
    consensus_input = ConsensusInput.parse_consensus_input(to_build)
    consensus_output = calculate_consensus(consensus_input=consensus_input)

    index_switches = SwitchData.parse_filtered_data("index", index_switches)
    major_switches = SwitchData.parse_filtered_data("major", major_switches)
    minor_switches = SwitchData.parse_filtered_data("minor", minor_switches)
    unique_switches = SwitchData.parse_filtered_data("unique", unique_switches)

    write_consensus(output_file_name,
                 consensus_input,
                 threshold,
                 consensus_output,
                 index_switches, major_switches, minor_switches, unique_switches)
    return consensus_output

def building_index(analyse,thold,dima_kmer):
    no_sup=[]
    to_build=[]
    low_sup=[]
    average=0
    highest=0

    for x in range(0,len(analyse.results)):
        average=average+analyse.results[x].support
        if analyse.results[x].support > highest:
            highest=analyse.results[x].support
        indexes=analyse.get_index('Index',x)
        if not indexes:
            indexes=x+1,'-'*dima_kmer,0,0
            no_sup.append(x+1)
            to_build.append(indexes)
        else:
            for k,v in indexes.items():

                a=v[0],k[0],analyse.results[x].support,v[1]
                to_build.append(a)
                if a[2] < thold:
                    low_sup.append(a[0])
    average = round(average / len(analyse.results))

    return to_build,no_sup,low_sup,highest,average

#below is to compare index and majors where they share an 8-mer

def unable_to_analyze(analyse,Filtered_Indexes,Filtered_Majors,dima_kmer,cons):
    unable_to_analyze=[]
    for x in range (0,len(analyse.results)-1):
        a = {'Index1':[],'Index2':[],'Major1':[],'Major2':[]}
        b = {}

        if x+2 in Filtered_Indexes or x+2 in Filtered_Majors:
            continue

        if analyse.results[x].support == 0 or analyse.results[x+1].support == 0:
            continue

        for i in analyse.results[x].variants:
            if i.motif_long == 'Index':
                a['Index1'].append(i.sequence[1::])
            if i.motif_long == 'Major':
                a['Major1'].append(i.sequence[1::])
        for i in analyse.results[x+1].variants:
            if i.motif_long == 'Index':
                a['Index2'].append(i.sequence[0:-1])
                b['Index']=i.sequence[-1]
            if i.motif_long == 'Major':
                a['Major2'].append(i.sequence[0:-1])
                b['Major']=i.sequence[-1]

        inter_i = set.intersection(set(a['Index1']), set(a['Major1']))
        inter_m = set.intersection(set(a['Index2']), set(a['Major2']))
        if inter_i == inter_m and inter_i:
            if cons[x+dima_kmer].isupper():
                if cons[x+dima_kmer] == b['Index']:
                    continue

                elif cons[x+dima_kmer] == b['Major']:
                    hits={}
                    fixa=analyse.get_index('Index',x)
                    fixb=analyse.get_index('Major',x+1)

                    pos1=list(fixa.values())[0][0]
                    seq1=list(fixa.keys())[0][0]
                    count1=list(fixa.values())[0][1]
                    tag1=list(fixa.keys())[0][1]

                    pos2=list(fixb.values())[0][0]
                    seq2=list(fixb.keys())[0][0]
                    count2=list(fixb.values())[0][1]
                    tag2=list(fixb.keys())[0][1]
                    hits.update({'({}){}'.format(seq2[:kmer],seq2[kmer:]):'seq with count {}'.format(count2),tag2:'tag at pos {}'.format(pos2)})
                    Filtered_Indexes.update({pos1+1:{'position': pos1, 'seq': '{}({}){}'.format(seq1[0],seq1[1:kmer+1],seq1[kmer+1:]),'count': count1,'tag':tag1, 'hits':hits}})
            else:
                unable_to_analyze.append('{} vs {}'.format(x+1,x+2))
    return unable_to_analyze

def check_which(which,Filtered_Indexes,Filtered_Majors,Filtered_Minors,Filtered_Uniques):
    if 'unique' not in which:
        Filtered_Uniques = {}
    if 'minor' not in which:
        Filtered_Minors = {}
    if 'major' not in which:
        Filtered_Majors = {}
    if 'index' not in which:
        Filtered_Indexes = {}

    return Filtered_Uniques,Filtered_Minors,Filtered_Majors,Filtered_Indexes
