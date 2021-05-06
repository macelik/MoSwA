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
    total_lines = len(consensus_input)
    for i in range(length_of_each_line + total_lines - 1):
        count = defaultdict(int)
        for j in range(max(0, i - length_of_each_line + 1), len(consensus_input)):
            if j > i:
                break
            alphabet = consensus_input[j].data[i - j]
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


def write_output(filepath,
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

    write_output(output_file_name,
                 consensus_input,
                 threshold,
                 consensus_output,
                 index_switches, major_switches, minor_switches, unique_switches)
    return consensus_output
