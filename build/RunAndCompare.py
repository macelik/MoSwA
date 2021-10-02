from dima import Dima
import os
import json,re

class GetMotifs():
    
    def __init__(self,fname,dima_kmer):
        print('Dima is running')
        self.results = Dima(sequences=fname, kmer_length=dima_kmer).run().results
        print('Dima is done')
    def get_index(self, motif, position):

        if motif == 'Minor':
            variants = self.results[position].get_minors('desc')
        else:
            variants = self.results[position].variants
        
        if not variants:
            return {}
        
        matches = {}
        
        for variant in variants:
            if variant.motif_long == motif:
                sequence = variant.sequence
                tag = '_'.join((motif, str(len(matches) + 1)))
                pos = position+1
                count = variant.count
                matches[sequence, tag] = pos,count

        if len(matches) == 1 and (motif == 'Index' or motif == 'Major'):
            matches[list(matches.keys())[0][0],motif] = list(matches.values())[0][0],list(matches.values())[0][1]
            s=list(matches.keys())[0][0],list(matches.keys())[0][1]
            matches.pop(s)
            

        return matches

class Compare():
    
    
    def __init__(self,tag,kmer):
        
        self.d = []
        self.key_value = dict()
        self.tag = tag
        self.kmer=kmer
        
    def search(self,analyse):
        for x in range(0, len(analyse.results)-1):
            self.key_value = dict()
            a,b=analyse.get_index(self.tag, x),analyse.get_index('Index', x+1)
            com2(a,b,self.d,self.key_value,self.kmer)

            a,b=analyse.get_index(self.tag, x),analyse.get_index('Major', x+1)
            com2(a,b,self.d,self.key_value,self.kmer)

            a,b=analyse.get_index(self.tag, x),analyse.get_index('Minor', x+1)
            com2(a,b,self.d,self.key_value,self.kmer)

            a,b=analyse.get_index(self.tag, x),analyse.get_index('Unique', x+1)
            com2(a,b,self.d,self.key_value,self.kmer)
        #return self.d
    
    def tag_motifs(self):
    
        for i in range(0,(len(self.d))):

            if self.d[i]['hits'] == {}:
                self.d[i]['tag']=self.d[i]['tag'] + ' x '
                #pos=pos+1
            elif not any(self.tag in key for key in self.d[i]['hits']):
                self.d[i]['tag']=self.d[i]['tag'] + ' x {}'.format(list(self.d[i]['hits'])[1])
        return self.d
    
    def filter_motifs(self):
        pos=0
        filtered={}
        
        for x in self.d:
            if x['position'] > pos:
                pos = x['position']

            if x['position'] == pos:

                if pos+1 in filtered:
                    pos = pos+1

                elif ' x ' in x['tag']:
                    filtered[x['position']+1] = x
        return filtered

def new_getem(data, m1):
    seq1=m1[0]
    pos1=data[m1][0]
    tag1=m1[1]
    count1=data[m1][1]
    return seq1,pos1,tag1,count1

def com2(a,b,d,key_value,kmer):
    for m1 in a:
        
        seq1, pos1, tag1, count1 = new_getem(a, m1)
        key = '{}-{}'.format(seq1[:], tag1)
        bla = dict()
        hits = dict()
        for m2 in b:
            if m1[0][1:kmer+1]==m2[0][:kmer]:
                seq2,pos2,tag2,count2=new_getem(b, m2)
                hits.update({'({}){}'.format(seq2[:kmer],seq2[kmer:]):'seq with count {}'.format(count2),tag2:'tag at pos {}'.format(pos2)})
        if key in key_value:
            data = d[key_value[key]]
            #print(key)
            #print(data)
            data['hits'].update(hits)
        else:
            bla.update({'position': pos1, 'seq': '{}({}){}'.format(seq1[0],seq1[1:kmer+1],seq1[kmer+1:]),'count': count1,'tag':tag1, 'hits':hits})
            d.append(bla)
            key_value[key] = len(d) - 1

def write_output(motif,raw_data,filtered_data):
    file_name='{}_raw.json'.format(motif)
    completeName=os.path.join(path, file_name)

    f = open(completeName.format(motif),'w')
    json.dump(raw_data,f)
    f.close()

    file_name='{}_filtered.json'.format(motif)
    completeName=os.path.join(path, file_name)

    f = open(completeName.format(motif),'w')
    json.dump(filtered_data,f)
    f.close()

def runit(which,dima_kmer,analyse):
    common_pos=[]
    lengthofswitches=dict()
    title={"Index":("Indexes","Raw_Indexes","Filtered_Indexes"),
              "Major":("Majors","Raw_Majors","Filtered_Majors"),
               "Minor":("Minors","Raw_Minors","Filtered_Minors"),
               "Unique":("Uniques","Raw_Uniques","Filtered_Uniques")
              }

    for x in title:

        if x.lower() in which:
            print("Analyzing",title[x][0])
            globals()[title[x][0]]=Compare(x,dima_kmer-1)
            globals()[title[x][0]].search(analyse)
            globals()[title[x][1]]=globals()[title[x][0]].tag_motifs()
            globals()[title[x][2]]=globals()[title[x][0]].filter_motifs()
            #write_output(x,globals()[title[x][1]],globals()[title[x][2]])
            lengthofswitches[x]=len(globals()[title[x][2]])
            common_pos.extend(globals()[title[x][2]].keys())
        else:
            globals()[title[x][0]]=Compare(x,dima_kmer-1)
            globals()[title[x][0]].search(analyse)
            globals()[title[x][1]]=globals()[title[x][0]].tag_motifs()
            globals()[title[x][2]]=globals()[title[x][0]].filter_motifs()
            lengthofswitches[x]="NA"
    return list(set(common_pos)),lengthofswitches,Raw_Indexes,Filtered_Indexes,Raw_Majors,Filtered_Majors,Raw_Minors,Filtered_Minors,Raw_Uniques,Filtered_Uniques

