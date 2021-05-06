def building_index(myd):
    to_build=[]
    no_sup=[]  ### added no_sup to show positions with no suppport
    for i in range(0,len(myd)):
        if myd[i]['variants'] == []:
            a=myd[i]['position'],'-'*len(myd[0]['variants'][0]['sequence']),myd[i]['supports'],0
            to_build.append(a)
            no_sup.append(myd[i]['position'])
        else:
            a=myd[i]['variants'][0]['position'],myd[i]['variants'][0]['sequence'],myd[i]['supports'],myd[i]['variants'][0]['count']
            to_build.append(a)
    return to_build,no_sup