
#import os,sys,inspect
#currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(currentdir)
#sys.path.insert(0,parentdir) 


def average_pos(myd):
	t_count=0
	e_pos=[]
	for i in range(len(myd)):
		if myd[i]['variants'] == []:
			e_pos.append(myd[i]['position'])
		for x in range(0,len(myd[i]['variants'])):
			t_count+=myd[i]['variants'][x]['count']

	average=t_count/len(myd)
	average="%.2f" % average
	return (average,e_pos)

def smthng():
	below=[]
	above=[]
	average_pos()
	for i in somethingelse:
		if indexes[i]['count'] > thold:
			above.append(i)
		else:
			below.append(i)
	return below,above





