import textwrap
def text1(myd,average,no_sup,low_sup,thold,common_pos):
    text1="""Your input alignment length is {} a.a with the highes supportt of 100 for a given position. The average number of sequences for a 

given kmer position is  {:.2f}. The number of positions that has no support (absence of kmer without indel gaps) is {} and the number of 

positions that has low support {} based on the threshold {} you submitted. There are {} number of unique switches out of {}"""
    
    text1=text1.format(len(myd),average,len(no_sup),len(low_sup),thold,len(set(common_pos)),(len(common_pos)))
    text1=textwrap.fill(text1, 136)
    text3="Below you are the positions that are without support\n\n {}\n\nAnd positions that has low support\n\n {}\n\n".format(textwrap.fill(', '.join(map(str,no_sup)), 152),textwrap.fill(', '.join(map(str,low_sup)), 152))
    return text1,text3

#def text3(which):
