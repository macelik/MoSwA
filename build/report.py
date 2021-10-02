import os

def write_report(average,
                thold,
                no_sup,
                low_sup,
                common_pos,
                lengthofswitches,
                unable_to_analyze,
                alignment_length,
                highest,
                path):
    summed = []
    if all(type(value) == int for value in lengthofswitches.values()):
        summed = sum(lengthofswitches.values())
    else:
        for x in lengthofswitches.values():
            if type(x) == int:
                summed.append(x)
        summed=sum(summed)


    replacements = {'alignment_length!':str(alignment_length), 'longest_support!':str(highest),
                'average!':str(average), 'len_no_sup!':str(len(no_sup)),'len_low_sup!':str(len(low_sup)),
                '!thold!':str(thold),'uniqueswitches!':str(len(common_pos)),'len_indexes!':str(lengthofswitches['Index']),
                'len_majors!':str(lengthofswitches['Major']),'len_minors!':str(lengthofswitches['Minor']),
                'len_uniques!':str(lengthofswitches['Unique']),'totalswitches!':str(summed),
                 '!UnableToAnalyze':'/'.join(map(str,unable_to_analyze))}

    with open('build/htmltemplate/index.html') as infile, open(os.path.join(path, "Report.html"), 'w') as outfile:
        for line in infile:
            for src, target in replacements.items():
                line = line.replace(src, target)
            outfile.write(line)
    return summed

def report_json(average,
                thold,
                no_sup,
                low_sup,
                common_pos,
                lengthofswitches,
                unable_to_analyze,
                alignment_length,
                highest,
                summed,
                dima_kmer,
                vis):

    vis['MoSwa_Output']['Report']['kmer_length']=dima_kmer
    vis['MoSwa_Output']['Report']['tresh_hold']=thold
    vis['MoSwa_Output']['Report']['alignment_length']=alignment_length
    vis['MoSwa_Output']['Report']['highest_support']=highest
    vis['MoSwa_Output']['Report']['average_support']=average
    vis['MoSwa_Output']['Report']['Positions_no_support']=len(no_sup)
    vis['MoSwa_Output']['Report']['Positions_low_support']=len(low_sup)
    vis['MoSwa_Output']['Report']['Total_switches']=lengthofswitches
    vis['MoSwa_Output']['Report']['Unique_switches']=len(common_pos)
    vis['MoSwa_Output']['Report']['Unable_to_analyze']='/'.join(map(str,unable_to_analyze))

    return vis



