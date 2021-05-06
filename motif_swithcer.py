#!/usr/bin/env python
from argparse import ArgumentParser
import os.path
import json
from build import ToBuild, ConsensusFromIndex, AnalyzeIndex, clustering, plotting, report, network_graph
from hunana import Hunana
rec_u=0
rec_m=0
bla={}
hits={}
d=[]
key_value = dict()  # seq-tag is stored as key, while index in the list is stored as value
indexes={}
majors={}
minors=[]
uniques=[]


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

parser = ArgumentParser(description="Motif finder")
parser.add_argument("-i", dest="filename", required=True,
					help="json file input", metavar="FILE",
					type=lambda x: is_valid_file(parser, x))
#parser.add_argument("-k", default=9, type=int, help="This is the k-mer size")
parser.add_argument("-m",
					choices=["index", "major", "minor", "unique","all"],
					nargs="+",
					default="all", type=str, help="motif sites to analyze")
parser.add_argument("-t", type=check_positive,required=True,
					help="Minimum number of seqeunces at a position")
args = parser.parse_args()
fname = args.filename
which= args.m
#kmer=args.k
thold=args.t

if fname.lower().endswith(('.fasta', '.fas', '.fa')):
	x = int(input("\nYour input is a fasta file, please provide k-mer length for Hunana: "))
	print('\nRunning Hunana')
	myd=Hunana(fname,kmer_len=x).run()

elif fname.lower().endswith(('.json')):
	with open(fname) as json_file:
		myd = json.load(json_file)

else:
	print('Your input has to be MSA or JSON file from Hunana output')
	sys.exit()


if 'all' in which:
	which='index','major','minor','unique'

#print (fname)

#print(myd)


def autoIncrementUniqeu():

	global rec_u
	pStart = 1 #adjust start value, if req'd
	pInterval = 1 #adjust interval value, if req'd
	if (rec_u == 0):
		rec_u = pStart
	else:
		rec_u = rec_u + pInterval
	return rec_u

def autoIncrementMinor():
	global rec_m

	pStart = 1 #adjust start value, if req'd
	pInterval = 1 #adjust interval value, if req'd
	if (rec_m == 0):
		rec_m = pStart
	else:
		rec_m = rec_m + pInterval
	return rec_m
#get rid of the unneeded info
def clean():
	for num in range (0,len(myd)):
		for i in range (len(myd[num]['variants'])):
			if "strain" in myd[num]['variants'][i]:
				myd[num]['variants'][i].pop("strain")
			if "motif_short" in myd[num]['variants'][i]:
				myd[num]['variants'][i].pop("motif_short")
			if "conservation" in myd[num]['variants'][i]:
				myd[num]['variants'][i].pop("conservation")
			if "id" in myd[num]['variants'][i]:
				myd[num]['variants'][i].pop("id")
			if "country" in myd[num]['variants'][i]:
				myd[num]['variants'][i].pop("country")
			if "host" in myd[num]['variants'][i]:
				myd[num]['variants'][i].pop("host")

#get index,major, minor and unique based on position
def get_index(pos):
	index={}
	for i in range (len(myd[pos]['variants'])):

		if myd[pos]['variants'][i]['motif_long']=='Index':

			index[myd[pos]['variants'][i]['sequence'],myd[pos]['variants'][i]['motif_long']]=myd[pos]['position'],myd[pos]['variants'][i]['count']

	return index

def get_major(pos):
	major={}
	for i in range (len(myd[pos]['variants'])):

		if myd[pos]['variants'][i]['motif_long']=='Major':

			major[myd[pos]['variants'][i]['sequence'],myd[pos]['variants'][i]['motif_long']]=myd[pos]['position'],myd[pos]['variants'][i]['count']

	return major



def get_minor(pos):
	minor={}
	global rec_m
	rec_m=0
	for i in range (len(myd[pos]['variants'])):

		if myd[pos]['variants'][i]['motif_long']=='Minor':

			minor[myd[pos]['variants'][i]['sequence'],myd[pos]['variants'][i]['motif_long']+'_'+ str(autoIncrementMinor())]=myd[pos]['position'],myd[pos]['variants'][i]['count']

	return minor

def get_unique(pos):
	unique={}
	global rec_u
	rec_u=0
	for i in range (len(myd[pos]['variants'])):

		if myd[pos]['variants'][i]['motif_long']=='Unique':

			unique[myd[pos]['variants'][i]['sequence'],myd[pos]['variants'][i]['motif_long']+'_'+ str(autoIncrementUniqeu())]=myd[pos]['position'],myd[pos]['variants'][i]['count']

	return unique

#getting the seq,count,tag separetly into variables
def new_getem(data, m1):
	seq1=m1[0]
	pos1=data[m1][0]
	tag1=m1[1]
	count1=data[m1][1]
	return seq1,pos1,tag1,count1

def com2(a,b):
	for m1 in a:

		seq1, pos1, tag1, count1 = new_getem(a, m1)
		key = '{}-{}'.format(seq1[:], tag1)
		#print(key)
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

def find_motif_index():
	for x in range (0,len(d)): #index working
		if d[x]['hits'] == {}:
			d[x]['tag']=d[x]['tag'] + ' x '
		elif 'Index' not in d[x]['hits'].keys():
			d[x]['tag']=d[x]['tag'] + ' x {}'.format(list(d[x]['hits'])[1])



def find_motif_major():
	for x in range (0,len(d)): #majors working
		if d[x]['hits'] == {}:
			d[x]['tag']=d[x]['tag'] + ' x '
		elif 'Major' not in d[x]['hits'].keys():
			d[x]['tag']=d[x]['tag'] + ' x {}'.format(list(d[x]['hits'])[1])


def find_motif_minor():
	pos=0  #working minors
	for i in range(0,(len(d))):

		if d[i]['position'] > pos:
			pos=d[i]['position']

		if d[i]['position'] == pos:

			if d[i]['hits'] == {}:
				d[i]['tag']=d[i]['tag'] + ' x '
				#pos=pos+1

			elif d[i]['tag'] not in d[i]['hits'].keys():
				if 'Minor' in list(d[i]['hits'].keys())[1]:
					continue
					#if int(list(d[i]['hits'].values())[0].split(" ")[-1]) != d[i]['count']:
						#d[i]['tag']=d[i]['tag'] + ' x {}'.format(list(d[i]['hits'])[1])
						#pos=pos+1

				else: #checking if count is not the same
					d[i]['tag']=d[i]['tag'] + ' x {}'.format(list(d[i]['hits'])[1])
					#pos=pos+1
				#pos=pos+1

def find_motif_unique():
	pos=0  #working unique
	for i in range(0,(len(d))):

		if d[i]['position'] > pos:
			pos=d[i]['position']

		if d[i]['position'] == pos:

			if d[i]['hits'] == {}:
				d[i]['tag']=d[i]['tag'] + ' x '


			elif d[i]['tag'] not in d[i]['hits'].keys():
				if 'Unique' in list(d[i]['hits'].keys())[1]:
					continue
				else:
					d[i]['tag']=d[i]['tag'] + ' x {}'.format(list(d[i]['hits'])[1])
				#pos=pos+1
			#elif 'Unique' not in str(d[i]['hits'].keys()):
				#d[i]['tag']=d[i]['tag'] + ' x {}'.format(list(d[i]['hits'])[1])
				#pos=pos+1



def write(i):
	i=str(i)
	f = open('{}.txt'.format(i),'w')
	json.dump(d,f)
	f.close()

def write_filtered(i):
	i=str(i)
	f = open('{}_filtered.txt'.format(i),'w')
	if i == 'index':
		json.dump(indexes,f)
		f.close()
	if i == 'major':
		json.dump(majors,f)
		f.close()
	if i == 'minor':
		global minor1
		minor1={}
		for item in minors:
			key=item['position']+1
			if key not in minor1:
				minor1[key] = [item]
			else:
				minor1[key].append(item)
		json.dump(minor1,f)
		f.close()

	if i == 'unique':
		global unique1
		unique1={}
		for item in uniques:
			key=item['position']+1
			if key not in unique1:
				unique1[key] = [item]
			else:
				unique1[key].append(item)
		json.dump(unique1,f)
		f.close()

def filter_indexes():
	for i in range(0,len(d)):
		if ' x ' in d[i]['tag']:
			indexes.update({d[i]['position']:d[i]})


def filter_majors():
	for i in range(0,len(d)):
		if ' x ' in d[i]['tag']:
			majors.update({d[i]['position']:d[i]})

def filter_minors():
	for i in range(0,len(d)):
		if ' x ' in d[i]['tag']:
			minors.append(d[i])

def filter_uniques():
	for i in range(0,len(d)):
		if ' x ' in d[i]['tag']:
			uniques.append(d[i])
			#uniques.update({d[i]['position']:d[i]})

def lets(which):
	global d,key_value,indexes,majors,minors,uniques,raw_index,raw_major,raw_minor,raw_unique

	for item in which:
		if item == 'index':
			d=list()
			key_value=dict()
			print ('scanning indexes')
			for i in range (0,len(myd)-1):
				a,b = get_index(i),get_index(i+1) #compare index with the following positions index
				com2(a,b)

				a,b = get_index(i),get_major(i+1) #compare index with the following positions major
				com2(a,b)

				a,b = get_index(i),get_minor(i+1) #compare index with the following positions minor
				com2(a,b)



				a,b = get_index(i),get_unique(i+1) #compare index with the following positions minor
				com2(a,b)
			raw_index=d
			find_motif_index()
			bla={}
			hits={}
			key_value = {}
			filter_indexes()
			repos_them={}
			for x in indexes:
				repos_them[x+1]=indexes[x]
			indexes=repos_them

			write(item)
			write_filtered(item)
		if item == 'major':
			d=list()
			key_value=dict()
			print ('scanning majors')
			for i in range (0,len(myd)-1):
				a,b = get_major(i),get_index(i+1) #compare index with the following positions index
				com2(a,b)

				a,b = get_major(i),get_major(i+1) #compare index with the following positions major
				com2(a,b)

				a,b = get_major(i),get_minor(i+1) #compare index with the following positions minor
				com2(a,b)



				a,b = get_major(i),get_unique(i+1) #compare index with the following positions minor
				com2(a,b)
			raw_major=d
			find_motif_major()
			bla={}
			hits={}
			key_value = {}
			filter_majors()
			repos_them={}
			for x in majors:
				repos_them[x+1]=majors[x]
			majors=repos_them

			write(item)
			write_filtered(item)


		if item == 'minor':
			print ('scanning minors')
			d=list()
			key_value=dict()
			for i in range (0,len(myd)-1):
				a,b = get_minor(i),get_index(i+1) #compare index with the following positions index
				com2(a,b)

				a,b = get_minor(i),get_major(i+1) #compare index with the following positions major
				com2(a,b)

				a,b = get_minor(i),get_minor(i+1) #compare index with the following positions minor
				com2(a,b)



				a,b = get_minor(i),get_unique(i+1) #compare index with the following positions minor
				com2(a,b)
			raw_minor=d
			find_motif_minor()
			bla={}
			hits={}
			key_value = {}
			filter_minors()
			write(item)
			write_filtered(item)
		if item == 'unique':
			print ('scanning unique')
			d=list()
			key_value=dict()
			for i in range (0,len(myd)-1):
				a,b = get_unique(i),get_index(i+1) #compare index with the following positions index
				com2(a,b)

				a,b = get_unique(i),get_major(i+1) #compare index with the following positions major
				com2(a,b)

				a,b = get_unique(i),get_minor(i+1) #compare index with the following positions minor
				com2(a,b)



				a,b = get_unique(i),get_unique(i+1) #compare index with the following positions minor
				com2(a,b)
			raw_unique=d
			find_motif_unique()
			bla={}
			hits={}
			key_value = {}
			filter_uniques()

	for item in which:
		if item == 'unique':
			write(item)
			write_filtered(item)


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
	dict1['Mi => +']=minorgain
	dict1['Mi => -']=minorloss
	dict1['U => +']=ugain
	dict1['U => -']=uniqueloss
	return dict1

def careful(myd,con):
	unable_to_analyze=[]
	for x in range (0,len(myd)-1):
		if list(get_index(x).keys()) == [] or list(get_index(x+1).keys()) == []:
			continue
		elif list(get_index(x).keys())[0][0][1::] == list(get_index(x+1).keys())[0][0][:-1]:
			if get_major(x+1) != {} and get_major(x) != {}:
				if list(get_major(x).keys())[0][0][1::] == list(get_major(x+1).keys())[0][0][:-1]:
					if list(get_index(x).keys())[0][0][1::] == list(get_major(x+1).keys())[0][0][:-1]:
						if con[x+9].isupper():
							if list(get_index(x+1).keys())[0][0][-1] != con[x+9]:
								hits={}
								pos1=list(get_index(x).values())[0][0]
								seq1=list(get_index(x).keys())[0][0]
								count1=list(get_index(x).values())[0][1]
								tag1=list(get_index(x).keys())[0][1]

								pos2=list(get_major(x+1).values())[0][0]
								seq2=list(get_major(x+1).keys())[0][0]
								count2=list(get_major(x+1).values())[0][1]
								tag2=list(get_major(x+1).keys())[0][1]
								hits.update({'({}){}'.format(seq2[:kmer],seq2[kmer:]):'seq with count {}'.format(count2),tag2:'tag at pos {}'.format(pos2)})
								indexes.update({pos1:{'position': pos1, 'seq': '{}({}){}'.format(seq1[0],seq1[1:kmer+1],seq1[kmer+1:]),'count': count1,'tag':tag1 + ' x ' + tag2, 'hits':hits}})
						else:
							unable_to_analyze.append('{} vs {}'.format(x+1,x+2))
	return indexes,unable_to_analyze

def analyze_index(indexes):
	indexloss=[]
	indextomajor=[]
	indextominor=[]
	indextounique=[]
	for i in indexes:
		if indexes[i]['hits'] == {}:
			indexloss.append(i)
			continue
		if indexes[i]['tag'].split(" x ",1)[1] == 'Major':
			indextomajor.append(i)
			continue
		if 'Minor' in indexes[i]['tag'].split(" x ",1)[1]:
			indextominor.append(i)
			continue
		if 'Unique' in indexes[i]['tag'].split(" x ",1)[1]:
			indextounique.append(i)
			continue

	return indexloss,indextomajor,indextominor,indextounique

def analyze_major(majors):
	majorloss=[]
	majortoindex=[]
	majortominor=[]
	majortounique=[]
	for i in majors:
		if majors[i]['hits'] == {}:
			majorloss.append(i)
			continue
		if majors[i]['tag'].split(" x ",1)[1] == 'Index':
			majortoindex.append(i)
			continue
		if 'Minor' in majors[i]['tag'].split(" x ",1)[1]:
			majortominor.append(i)
			continue
		if 'Unique' in majors[i]['tag'].split(" x ",1)[1]:
			majortounique.append(i)
			continue
	return majorloss,majortoindex,majortominor,majortounique

def analyze_minor(minors):
	minorloss=[]
	minortoindex=[]
	minortomajor=[]
	minortominor=[]
	minortounique=[]
	for i in minors:
		if i['hits'] == {}:
			minorloss.append(i['position']+1)
			continue
		if i['tag'].split(" x ",1)[1] == 'Index':
			minortoindex.append(i['position']+1)
			continue
		if i['tag'].split(" x ",1)[1] == 'Major':
			minortomajor.append(i['position']+1)
			continue
		if 'Minor' in i['tag'].split(" x ",1)[1]:
			minortominor.append(i['position']+1)
			continue
		if 'Unique' in i['tag'].split(" x ",1)[1]:
			minortounique.append(i['position']+1)
			continue
	return minorloss,minortoindex,minortomajor,minortominor,minortounique

def analyze_unique(uniques):
	uniqueloss=[]
	uniquetoindex=[]
	uniquetomajor=[]
	uniquetominor=[]
	uniquetounique=[]
	for i in uniques:
		if i['hits'] == {}:
			uniqueloss.append(i['position']+1)
			continue
		if i['tag'].split(" x ",1)[1] == 'Index':
			uniquetoindex.append(i['position']+1)
			continue
		if i['tag'].split(" x ",1)[1] == 'Major':
			uniquetomajor.append(i['position']+1)
			continue
		if 'Minor' in i['tag'].split(" x ",1)[1]:
			uniquetominor.append(i['position']+1)
			continue
		if i['tag'].split(" x ",1)[1] == 'Unique':
			uniquetounique.append(i['position']+1)
			continue
	return uniqueloss,uniquetoindex,uniquetomajor,uniquetominor,uniquetounique


def i_splits_mergers(indexes,majortoindex,minortoindex,uniquetoindex):
	mtoi,mitoi,utoi,igain=[],[],[],[]
	for i in indexes:
		if i in majortoindex:
			mtoi.append(i)
			continue
		elif i in minortoindex:
			mitoi.append(i)
		if list(get_index(i-1).keys()) == []:
			continue
		elif list(get_index(i-1).keys())[0][0][:-1] in str(get_minor(i-2).keys()):
			mitoi.append(i)
		elif i in uniquetoindex:
			utoi.append(i)
		elif i not in mtoi + mitoi + utoi:
			igain.append(i)
	merge_mtoi,merge_mitoi,merge_utoi=[],[],[]
	for x in majortoindex:
		if x not in mtoi:
			merge_mtoi.append(x)
	for x in set(minortoindex):
		if x not in mitoi:
			merge_mitoi.append(x)
	for x in uniquetoindex:
		if x not in utoi:
			merge_utoi.append(x)
	return mtoi,mitoi,utoi,igain,merge_mtoi,merge_mitoi,merge_utoi

def m_splits_mergers(majors,indextomajor,minortomajor,uniquetomajor):

	split_itom,itom,mitom,merge_mitom,utom,merge_utom,mgain=([] for i in range(7))

	for i in majors:
		if i in indextomajor:
			itom.append(i)
			continue
		if get_major(i-1) != {}:
			if list(get_index(i-2).keys())[0][0][1::] == list(get_major(i-1).keys())[0][0][:-1]:
				split_itom.append(i)
				continue
		if i in minortomajor:
			mitom.append(i)
			continue
		if len(get_major(i-1)) != 0:
			if list(get_major(i-1).keys())[0][0][:-1] in str(get_unique(i-2).keys()):
				mitom.append(i)
				continue
		if i in uniquetomajor:
			utom.append(i)
			continue
		if i not in itom+mitom+utom+split_itom:
			mgain.append(i)
			continue
	for x in minortomajor:
		if x not in mitom:
			merge_mitom.append(x)
	for x in uniquetomajor:
		if x not in utom:
			merge_utom.append(x)

	return split_itom,itom,mitom,merge_mitom,utom,merge_utom,mgain


def mi_splits_mergers(minors,indextominor,majortominor,uniquetominor):

	split_itomi,itomi,split_mtomi,mtomi,utomi,merge_utomi,migain,minortominor=([] for i in range(8))
	count=0
	for i in minors:

		pos=i['position']
		tag=i['tag'].split(" x ",1)[0]
		seq2=''
		for x in list(get_minor(pos).keys()):
			if x[1] == tag:
				seq2=x[0][:-1]

		if i['position']+1 in indextominor:
			#itomi.append(i['position'])
			continue

		elif seq2 != '' and seq2 == list(get_index(pos-1).keys())[0][0][1::]:
			split_itomi.append(pos+1)
			continue

		elif i['position']+1 in majortominor:
			#mtomi.append(i['position'])
			continue

		elif seq2 != '' and seq2 == list(get_major(pos-1).keys())[0][0][1::]:
			split_mtomi.append(i['position']+1)
			continue

		elif seq2 != '':
			for y in list(get_minor(pos-1).keys()):
				if seq2 == y[0][1::]:
					minortominor.append(pos+1)
					break
		elif i['position']+1 in uniquetominor:
			utomi.append(i['position']+1)
			#print(i['position'],i['tag'],seq2)
			continue

	#Minor Gain
	count=0
	minorgain=[]
	for x in range (0,len(myd)-1):
		#if len(get_minor(x+1)) < len(get_minor(x)):
			#continue
		#else:
		ai,bi,ci = {},{},{}
		a=get_minor(x)
		for keys in a:
			ai[keys[1]]=keys[0][1::]
		b=get_minor(x+1)
		for keys in b:
			bi[keys[1]]=keys[0][:-1]
		c=get_unique(x)
		for keys in c:
			ci[keys[1]]=keys[0][1::]

		for v in bi.values():
			if v in ai.values():
				continue
			if list(get_index(x).keys()) ==[]:
				continue
			if v == list(get_index(x).keys())[0][0][1::]:
				#print('from index',v,x)
				continue
			if list(get_major(x).keys()) == []:
				continue
			if v == list(get_major(x).keys())[0][0][1::]:
				continue
			if v in ci.values():
				continue
			else:
				minorgain.append(x+2)
				break
	itomi=list(set(itomi))
	split_itomi=list(set(split_itomi))
	split_mtomi=list(set(split_mtomi))
	mtomi=list(set(mtomi))
	split_mtomi=list(set(split_mtomi))
	utomi=list(set(utomi))
	minortominor=list(set(minortominor))
	for x in uniquetominor:
		if x not in utomi:
			merge_utomi.append(x)

	return split_itomi,itomi,split_mtomi,mtomi,utomi,merge_utomi,migain,minortominor,minorgain


def u_splits_mergers(uniques,indextounique,majortounique,minortounique):

	split_itou,itou,split_mtou,mtou,split_mitou,mitou,ugain,utou=([] for i in range(8))
	count=0
	for i in uniques:

		pos=i['position']
		tag=i['tag'].split(" x ",1)[0]
		seq2=''
		for x in list(get_unique(pos).keys()):
			if x[1] == tag:
				seq2=x[0][:-1]

		if i['position']+1 in indextounique:
			#itomi.append(i['position'])
			continue

		elif seq2 != '' and seq2 == list(get_index(pos-1).keys())[0][0][1::]:
			split_itou.append(pos+1)
			continue

		elif i['position']+1 in majortounique:
			#mtomi.append(i['position'])
			continue

		elif seq2 != '' and seq2 == list(get_major(pos-1).keys())[0][0][1::]:
			split_mtou.append(i['position']+1)
			continue


		elif i['position']+1 in minortounique:
			#utomi.append(i['position'])
			#print(i['position'],i['tag'],seq2)
			continue

		elif seq2 != '':
			for y in list(get_minor(pos-1).keys()):
				if seq2 == y[0][1::]:
					split_mitou.append(pos+1)
					break
		if tag != '' and pos+1 not in (indextounique+majortounique+minortounique+split_itou+split_mtou+split_mitou):
			#print('here')
			for a in list(get_unique(pos-1).keys()):
				if seq2 == a[0][1::]:
					utou.append(i['position']+1)
					break

			if seq2 != '' and pos+1 not in utou:
				ugain.append(pos+1)
				continue
		#if seq2 == '':
			#print(tag,pos)
			#count += 1


	split_itou=list(set(split_itou))
	split_mtou=list(set(split_mtou))
	split_mitou=list(set(split_mitou))
	utou=list(set(utou))
	minortounique=list(set(minortounique))

	count=0
	uniquegain=[]
	for x in range (0,len(myd)-1):
		#if len(get_minor(x+1)) < len(get_minor(x)):
			#continue
		#else:
		ai,bi,ci = {},{},{}
		a=get_unique(x)
		for keys in a:
			ai[keys[1]]=keys[0][1::]
		b=get_unique(x+1)
		for keys in b:
			bi[keys[1]]=keys[0][:-1]
		c=get_minor(x)
		for keys in c:
			ci[keys[1]]=keys[0][1::]

		for v in bi.values():
			if v in ai.values():
				continue
			if list(get_index(x).keys()) == []:
				continue
			if v == list(get_index(x).keys())[0][0][1::]:
				#print('from index',v,x)
				continue
			if list(get_major(x).keys()) == []:
				continue
			if v == list(get_major(x).keys())[0][0][1::]:
				continue
			if v in ci.values():
				continue
			else:
				uniquegain.append(x+1)
				break

	return split_itou,itou,split_mtou,mtou,split_mitou,mitou,ugain,utou,uniquegain

def align_i():
	import csv
	from Bio import pairwise2
	from Bio.Align.substitution_matrices import load as ld
	matrix = ld(("PAM30"))
	gap_open = -10
	gap_extend = -0.5
	align={}

	for x in indextomajor:
		align[x]={'Index':list(get_index(x-1).keys())[0][0],'Major':list(get_major(x-1).keys())[0][0]}
		a=align[x]['Index']
		b=align[x]['Major']
		alns = pairwise2.align.globalds(a, b, matrix, gap_open, gap_extend)
		align[x]['score']=alns[0][2]

	for x in indextominor:
		pos=x
		tag=list(indexes[x].values())[3].split(" ")[-1]
		seq1=list(get_index(x-1).keys())[0][0]
		for i in list(get_minor(x-1).keys()):
			if i[1] == tag:
				seq2 = i[0]
				align[x]={'Index':seq1,str(tag):seq2}
				a=align[x]['Index']
				b=align[x][str(tag)]

				alns = pairwise2.align.globalds(a, b, matrix, gap_open, gap_extend)
				align[x]['score']=alns[0][2]

	for x in indextounique:
		pos=x
		tag=list(indexes[x].values())[3].split(" ")[-1]
		seq1=list(get_index(x-1).keys())[0][0]
		for i in list(get_unique(x-1).keys()):
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

	with open("Index_Switches_Alignment.csv", "w") as csv_file:
		csvwriter = csv.writer(csv_file)
		csvwriter.writerow(['Position', 'Index', 'Motif', 'Score'])

		for x in fix_align:
			if x in majortoindex:
				tag='Major'
			elif x in minortoindex:
				tag='Minor'
			elif x in uniquetoindex:
				tag='Uniqe'
			csvwriter.writerow([x, fix_align[x]['Index']+' ({})'.format(tag),fix_align[x][list(fix_align[x])[1]]+' ({})'.format(list(fix_align[x])[1]), fix_align[x]['score']])

	return align,fix_align

kmer=len(myd[0]['variants'][0]['sequence'])-1
clean()
print("Selected motif search sites: {} \nMotif search kmer ({}-1): {}".format(args.m,kmer+1,kmer))
print("\nMoSwA is running..")
lets(which)
common_pos=[]

for item in which:
	if item == 'index':
		print ('\nThe number of switches found in index is {}'.format(len(indexes)))

		for i in indexes:
			common_pos.append(i)
	if item == 'major':
		print ('The number of switches found in major is {}'.format(len(majors)))

		for i in majors:
			common_pos.append(i)
	if item == 'minor':
		print ('The number of switches found in minor is {}'.format(len(minors)))

		for i in minors:
			common_pos.append(i['position'])
	if item == 'unique':
		print ('The number of switches found in unique is {}'.format(len(uniques)))

		for i in uniques:
			common_pos.append(i['position'])

print ('Unique number of switches: {} '.format(len(list(set(common_pos)))))
#minor1=minors
#unique1=uniques
to_build,no_sup=ToBuild.building_index(myd)
print("\nFollowing positions have no support thus not included in the analysis {}".format(', '.join(map(str,no_sup))))

indexloss,indextomajor,indextominor,indextounique=analyze_index(indexes)
cons = ConsensusFromIndex.ConsensusInput.parse_consensus_input(to_build)
cons = consensus_output = ConsensusFromIndex.calculate_consensus(consensus_input=cons)
indexes,unable_to_analyze=careful(myd,cons)
write_filtered('index')
majorloss,majortoindex,majortominor,majortounique=analyze_major(majors)
minorloss,minortoindex,minortomajor,minortominor,minortounique=analyze_minor(minors)
uniqueloss,uniquetoindex,uniquetomajor,uniquetominor,uniquetounique=analyze_unique(uniques)
mtoi,mitoi,utoi,igain,merge_mtoi,merge_mitoi,merge_utoi=i_splits_mergers(indexes,majortoindex,minortoindex,uniquetoindex)
split_itom,itom,mitom,merge_mitom,utom,merge_utom,mgain=m_splits_mergers(majors,indextomajor,minortomajor,uniquetomajor)
split_itomi,itomi,split_mtomi,mtomi,utomi,merge_utomi,migain,minortominor,minorgain=mi_splits_mergers(minors,indextominor,majortominor,uniquetominor)
split_itou,itou,split_mtou,mtou,split_mitou,mitou,ugain,utou,uniquegain=u_splits_mergers(uniques,indextounique,majortounique,minortounique)
if 'unique' not in which:
	unique1 = {}
if 'minor' not in which:
	minor1 = {}
if 'major' not in which:
	major = {}
if 'index' not in which:
	indexes = {}
ConsensusFromIndex.build_consensus(indexes, majors, minor1, unique1, to_build, thold, "output.fasta")

average,e_pos=AnalyzeIndex.average_pos(myd)


if 'unique' not in which:
	uniquegain,uniqueloss,minor1 = [],[],[]
	if 'minor' not in which:
		minorgain,minorloss,uniquetominor,minortounique = [],[],[],[]
	if 'major' not in which:
		majorgain,majorloss,uniquetomajor,majortounique = [],[],[],[]
	if 'index' not in which:
		indexgain,indexloss,uniquetoindex,indextounique = [],[],[],[]

if 'minor' not in which:
	minorgain,minorloss = [],[]
	if 'unique' not in which:
		uniquegain,uniqueloss,uniquetominor,minortounique = [],[],[],[]
	if 'major' not in which:
		majorgain,majorloss,minortomajor,majortominor= [],[],[],[]
	if 'index' not in which:
		indexgain,indexloss,indextominor,minortoindex = [],[],[],[]
if 'major' not in which:
	majorgain,majorloss = [],[]
	if 'unique' not in which:
		uniquegain,uniqueloss,majortounique,uniquetomajor = [],[],[],[]
	if 'minor' not in which:
		minorloss,minorgain,majortominor,minortomajor = [],[],[],[]
	if 'index' not in which:
		indexgain,indexloss,indextomajor,majortoindex = [],[],[],[]
if 'index' not in which:
	indexgain,indexloss = [],[]
	if 'unique' not in which:
		uniquegain,uniqueloss,indextounique,uniquetoindex = [],[],[],[]
	if 'minor' not in which:
		minorloss,minorgain,indextominor,minortoindex = [],[],[],[]
	if 'major' not in which:
		majorloss,majorgain,indextomajor,majortoindex = [],[],[],[]




#print("\nnumber of lost positions",len(indexlost))
dict1=create_dic()
print('\nClustering...')
low_no_sup,no_sup,low_sup=clustering.low_no_sup(myd,thold)
topla=clustering.uniqpos_no_low_support(common_pos,low_no_sup)
cluster,r_cluster=clustering.cluster_it(topla)
clustering.CnandCs(r_cluster,kmer)

hots,r_hots=clustering.hotspots(r_cluster)
plotin=plotting.plotit(r_cluster,kmer)

fig=plotting.draw(r_hots,plotin)
report.write_report(myd,average,no_sup,low_sup,thold,topla,indexes,majors,minors,uniques,indextomajor,majortoindex,indextominor,minortoindex,
				 indextounique,uniquetoindex,igain,indexloss,majortominor,minortomajor,majortounique,uniquetomajor,mgain,majorloss,minortounique,
				 uniquetominor,minorgain,minorloss,uniquegain,uniqueloss,common_pos,unable_to_analyze)
#text1,text3=texts.text1(myd,average,no_sup,low_sup,thold,common_pos)
#fig1=plotting.firstpage(text1,general_info,text3)
#plotting.write_pdf(fig,fig1)
sandm={'split_itom':split_itom,'split_itomi':split_itomi,'split_itou':split_itou,
	   'split_mtomi':split_mtomi,'split_mtou':split_mtou,'split_mitou':split_mitou,
	   'merge_mtoi':merge_mtoi,'merge_mitoi':merge_mitoi,'merge_utoi':merge_utoi,
	   'merge_mitom':merge_mitom,'merge_utom':merge_utom,'merge_utomi':merge_utomi}

print('\nPlotting network graph of switched positions...')
repos,D,connection=network_graph.network(dict1)
s_edges,s_nodes = network_graph.plotly_edges(repos,D,connection,sandm)
data=network_graph.plotly_info(s_edges,s_nodes)
network_graph.out_plotly(data)

fin = open("Motif_Plots.html", "rt")
data = fin.read()
data = data.replace("g.tx=\"Aa\"","g.tx=\"\"")
fin.close()
fin = open("Motif_Plots.html", "wt")
fin.write(data)
fin.close()
#align=align_i()
#print(align)
align,fix_align=align_i()
