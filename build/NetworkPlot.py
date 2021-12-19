import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout
import random,os
#from graphviz import Digraph
import plotly.graph_objs as go
import plotly.express as px
from plotly.offline import download_plotlyjs, init_notebook_mode, iplot

fr={'I => M': 'M => I', 'I => Mi': 'Mi => I', 'I => U': 'U => I', 'I => -': 'I => +', 'M => Mi': 'Mi => M', 'M => U': 'U => M', 'M => -': 'M => +', 'Mi => U': 'U => Mi', 'Mi => -': 'Mi => +', 'U => -': 'U => +'}
legend_d={'I => M': 'g1', 'M => I': 'g2', 'I => Mi': 'g3', 'Mi => I': 'g4', 'I => U': 'g5', 'U => I': 'g6', 'I => +': 'g7', 'I => -': 'g8', 'M => Mi': 'g9', 'Mi => M': 'g10', 'M => U': 'g11', 'U => M': 'g12', 'M => +': 'g13', 'M => -': 'g14', 'Mi => U': 'g15', 'U => Mi': 'g16', 'Mi => +': 'g17', 'Mi => -': 'g18', 'U => +': 'g19', 'U => -': 'g20'}

def check_arg(dict1,which):

    if 'unique' not in which:
        dict1['U => I']=[]
        dict1['U => M']=[]
        dict1['U => Mi']=[]
        dict1['U => +']=[]
        dict1['U => -']=[]
    if 'minor' not in which:
        dict1['Mi => I']=[]
        dict1['Mi => M']=[]
        dict1['Mi => U']=[]
        dict1['Mi => +']=[]
        dict1['Mi => -']=[]
    if 'major' not in which:
        dict1['M => I']=[]
        dict1['M => Mi']=[]
        dict1['M => U']=[]
        dict1['M => +']=[]
        dict1['M => -']=[]
    if 'index' not in which:
        dict1['I => M']=[]
        dict1['I => Mi']=[]
        dict1['I => U']=[]
        dict1['I => +']=[]
        dict1['I => -']=[]
    return dict1

def network(dict1):
    dict2={k: v for k, v in dict1.items() if v != []}
    D=nx.from_dict_of_lists(dict2)

    connection={}
    for node,adjacencies in enumerate(D.adjacency()):
        if type(adjacencies[0]) == int:
            connection[adjacencies[0]]=list(adjacencies[1].keys())

    color_map=[]
    for node in D:
        if type(node) == int:
            color_map.append('lightblue')
        else:
            color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
            color_map.append(color)

    nodeSize = 2500
    nodeColorList = color_map

    prog = 'dot'
    args = '-Gnodesep=50 -Grank=source -Gpad=0.5 -Gordering=out'
    root = None
    pos  = graphviz_layout(D, prog = prog, root = root, args = args)

    repos_PN={}
    repos_CN={}

    for key in pos:
        if type(key) != int:
            repos_PN[key]=pos[key]
        else:
            repos_CN[key]=pos[key]
    PN_sorted=dict(sorted(repos_PN.items(), key=lambda item: item[1][0]))
    CN_sorted=dict(sorted(repos_CN.items(), key=lambda item: item[1][0]))

    start=list(PN_sorted.values())[0][0]
    end=list(PN_sorted.values())[-1][0]
    space=(end-start)/((len(PN_sorted)/2))


    for k1,k2 in zip(repos_PN,PN_sorted):
        repos_PN[k1]=PN_sorted[k2]
    for k1,k2 in zip(repos_CN,CN_sorted):
        repos_CN[k1]=CN_sorted[k2]

    for k,v in fr.items():
        if k in PN_sorted:
            repos_PN[k]=(start,90)
        if v in PN_sorted:
            repos_PN[v]=(start,-90)
        if k in repos_PN or v in repos_PN:
            start=start+space

    repos={**repos_PN, **repos_CN}

    return repos,D,connection

def plotly_edges(repos,D,connection,sandm):

    colorscale2 = px.colors.qualitative.Vivid

    for n, p in repos.items():
        D.nodes[n]['pos'] = p

    s_edges=[]
    for edge in D.edges():
        edge_trace = go.Scatter(
        x=[],
        y=[],
        line=dict(width=0.5,color='#888'),
        hoverinfo='none',
        mode='lines')
        x0, y0 = D.nodes[edge[0]]['pos']
        x1, y1 = D.nodes[edge[1]]['pos']
        edge_trace['x'] = tuple([x0, x1, None])
        edge_trace['y'] = tuple([y0, y1, None])
        edge_trace['name']=edge[0]
        edge_trace['legendgroup']=legend_d[edge[0]]
        edge_trace['showlegend']=False
        #edge_trace['hovertext']=edge[0] + ' ' + str(edge[1])
        #edge_trace['line']['text']= edge[0] + ' ' + str(edge[1])
        if edge[0] in list(fr.keys()):
            if edge[0] in ['I => +','M => +','Mi => +','U => +']:
                #print('this')
                s_edges.append(edge_trace)
                #continue
            elif edge[0] == 'I => M':
                    #print('itom')
                if edge[1] in sandm['split_itom']:
                    #print('SPLIT')
                    edge_trace['line']['color'] = 'lightblue'
                    #continue
                else:
                    s_edges.append(edge_trace)
                    #continue

            elif edge[0] == 'I => Mi':
                if edge[1] in sandm['split_itomi']:
                    edge_trace['line']['color'] = 'lightblue'
                else:
                    s_edges.append(edge_trace)
            elif edge[0] == 'I => U':
                if edge[1] in sandm['split_itou']:
                    edge_trace['line']['color'] = 'lightblue'
                else:
                    s_edges.append(edge_trace)
            elif edge[0] == 'M => Mi':
                if edge[1] in sandm['split_mtomi']:
                    edge_trace['line']['color'] = 'lightblue'
                else:
                    s_edges.append(edge_trace)
            elif edge[0] == 'M => U':
                if edge[1] in sandm['split_mtou']:
                    edge_trace['line']['color'] = 'lightblue'
                else:
                    s_edges.append(edge_trace)
            elif edge[0] == 'Mi => U':
                if edge[1] in sandm['split_mitou']:
                    edge_trace['line']['color'] = 'lightblue'
                else:
                    s_edges.append(edge_trace)

        elif edge[0] in list(fr.values()):
            if edge[0] in ['I => -','M => -','Mi => -','U => -']:
                s_edges.append(edge_trace)

            elif edge[0] == 'M => I':
                if edge[1] in sandm['merge_mtoi']:
                    edge_trace['line']['color'] = 'coral'
                else:
                    s_edges.append(edge_trace)
            elif edge[0] == 'Mi => I':
                if edge[1] in sandm['merge_mitoi']:
                    edge_trace['line']['color'] = 'coral'
                else:
                    s_edges.append(edge_trace)

            elif edge[0] == 'U => I':
                if edge[1] in sandm['merge_utoi']:
                    edge_trace['line']['color'] = 'coral'
                else:
                    s_edges.append(edge_trace)

            elif edge[0] == 'Mi => M':
                if edge[1] in sandm['merge_mitom']:
                    edge_trace['line']['color'] = 'coral'
                else:
                    s_edges.append(edge_trace)

            elif edge[0] == 'U => M':
                if edge[1] in sandm['merge_utom']:
                    edge_trace['line']['color'] = 'coral'
                else:
                    s_edges.append(edge_trace)

            elif edge[0] == 'U => Mi':
                if edge[1] in sandm['merge_utomi']:
                    edge_trace['line']['color'] = 'coral'
                else:
                    s_edges.append(edge_trace)


            #elif edge[1] in (merge_mtoi+merge_mitoi+merge_utoi+merge_mitom+merge_utom+merge_utomi):
                #edge_trace['line']['color'] = 'lightgray'
        s_edges.append(edge_trace)

    num={1:'\u00b9',2:'\u00b2',3:'\u00b3',4:'\u2074',5:'\u2075',6:'\u2076',7:'\u2077',8:'\u2078',9:'\u2079',10:'\u00b9\u2070',11:'\u00b9\u00b9',12:'\u00b9\u00b2',13:'\u00b9\u00b3',14:'\u00b9\u2074',15:'\u00b9\u2075',16:'\u00b9\u2076',17:'\u00b9\u2077',18:'\u00b9\u2078',19:'\u00b9\u2079',20:'\u00b2\u2070'}
    s_nodes=[]
    count=0
    m=0
    n=10
    for _, adjacencies in enumerate(D.adjacency()):
        if m < len(adjacencies[1]):
            m = len(adjacencies[1])
        if n > len(adjacencies[1]):
            n = len(adjacencies[1])
    for _, adjacencies in enumerate(D.adjacency()):
        if type(adjacencies[0]) != int:
            con=''
            node_trace = go.Scatter(
            x=[],
            y=[],
            text=[],
            hovertext=[],
            marker=dict(
                #showscale=True,
                colorscale=colorscale2,
                cmin=n,
                cmax=m,
                size=15,
                color=[],
                line=dict(width=0)),
            mode='markers+text',
            hoverinfo='text')


            count += 1
            x, y = D.nodes[adjacencies[0]]['pos']
            node_trace['x'] = tuple([x])
            node_trace['y'] = tuple([y])
            node_trace['name'] = num[count] + adjacencies[0]
            node_trace['legendgroup']=legend_d[adjacencies[0]]
            node_trace['showlegend']=True

            node_info = str(adjacencies[0]) +' number of switches: '+str(len(adjacencies[1]))
            node_trace['hovertext']+=tuple([node_info])
            node_trace['text']=str(count)
            node_trace['marker']['color']=[len(adjacencies[1])]
            s_nodes.append(node_trace)
            for a in list(adjacencies[1].keys()):
                con=''
                node_trace = go.Scatter(
                x=[],
                y=[],
                text=[],
                hovertext=[],
                marker=dict(
                    size=15,
                    colorscale='Cividis',
                    cmin=1,
                    cmax=20,
                    color=[]),
                mode='markers+text',
                hoverinfo='text')
                x, y = D.nodes[a]['pos']
                node_trace['x'] = tuple([x])
                node_trace['y'] = tuple([y])
                node_trace['name'] = adjacencies[0]
                node_trace['legendgroup']=legend_d[adjacencies[0]]
                node_trace['showlegend']=False
                node_trace['marker']['color']=len(connection[a]),
                if len(connection[a]) > 1:

                    con= ', '.join(connection[a])
                    node_info = str(a) +' Multiple Switches b/w: '+ con
                    node_trace['hovertext']+=tuple([node_info])
                else:
                    con= ', '.join(connection[a])
                    node_info = str(a) +' Single switch b/w: '+ con
                    node_trace['hovertext']+=tuple([node_info])
                s_nodes.append(node_trace)


    return s_edges,s_nodes

def plotly_info(s_edges,s_nodes):
    data=[]
    for trace in s_edges:
        data.append(trace)
    for trace in s_nodes:
        data.append(trace)
    return data

def out_plotly(datas,path):
    fig = go.FigureWidget(data=datas,
             layout=go.Layout(
                #autosize=False,
                title='<br>Motif Switch Connections',
                titlefont=dict(size=16),
                showlegend=True,
                hovermode='closest',
                #clickmode='lasso+select',
                margin=dict(b=20,l=5,r=5,t=40),

                annotations=[ dict(
                    text="No. of connections",
                    showarrow=False,
                    #xref="paper", yref="paper",
                    xanchor='left',y=-110) ],
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))

    fig.write_html(os.path.join(path, "Motif_Plots.html"))

def fix_html(path):
    fin = open(os.path.join(path, "Motif_Plots.html"), "rt")
    data = fin.read()
    data = data.replace("g.tx=\"Aa\"","g.tx=\"\"")
    fin.close()
    fin = open(os.path.join(path, "Motif_Plots.html"), "wt")
    fin.write(data)
    fin.close()
