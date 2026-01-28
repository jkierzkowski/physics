import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.animation import FuncAnimation
import subprocess
import time
import ast

def adjacency_matrix(all_nodes, edges):
  A = np.zeros((len(edges), len(all_nodes)))
  edge_indices = np.array([[e[0], e[1]] for e in edges])
  A[np.arange(len(edges))[:, None], edge_indices] = 1
  return A

def length(n1, n2):
    return np.sqrt((n1[0] - n2[0]) ** 2 + (n1[1] - n2[1]) ** 2)

def change_neighbors(neighbors,nodes):
    return [[nodes.index(n) for n in n_e] for n_e in neighbors if len(n_e) > 1]


def lattice(n = 5,m = 5):

    triag = nx.triangular_lattice_graph(m,n,with_positions=True)
    nodes = list(triag.nodes())
    edges = [(nodes.index(e0),nodes.index(e1)) for e0,e1 in list(triag.edges())]
    pos = [[it,p[1]] for it,p in enumerate(list(triag.nodes.data("pos")))]
    neighbors = [list(triag.neighbors(node)) for node in triag.nodes()]
 
    return edges,nodes,pos,change_neighbors(neighbors,nodes)

def merging_matrix(D,pos_dict,edges):
    out,idx = [],[]
    for i, d0 in enumerate(D[:-1]):
        for j, d1 in enumerate(D[i+1:]):
            e2 = np.argwhere((np.sign(d0)-np.sign(d1))!=0)
            if np.count_nonzero(np.sign(d0)-np.sign(d1))<3 and (e2[0],e2[1]) in edges:
                l = length(pos_dict[e2[0,0]],pos_dict[e2[1,0]])
                m = d0[np.nonzero(d0)][0]+d1[np.nonzero(d1)][0]-l
                idx.append((i,j+i+1))
                out.append(m)
    return out, idx 
def save_to_file(D,path):
    f = open(path,"w")
    for d in D:
        if type(d)==list:
            for x in d:
                f.write(str(x).replace("(", "").replace(")", "")+' ')
            f.write('\n')
        else:
            f.write(str(d).replace(",","").replace("(", "").replace(")", ""))
            f.write('\n')
    f.close()
def import_edges(path="infinite_edges.txt"):
    deleted_edges=[]
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in  lines:
            edge_data = ast.literal_eval(line)
            deleted_edges.append((edge_data[0],edge_data[1]))
    return deleted_edges

def import_timeline(path="timeline.txt"):
    graphs = []
    # Open the file
    with open(path, 'r') as file:
        lines = file.readlines()

    G = nx.Graph()
    for line in lines:
        line = line.strip()
        if not line.startswith('-'):
            edge_data = ast.literal_eval(line)
            G.add_edge(edge_data[0][0], edge_data[0][1], diameter=edge_data[1])
        elif len(G.edges)>0:
            graphs.append(G)
            G = nx.Graph()
    return graphs


def gif_gen(fp, latt_list, pos,inf_edges,merged_edges):
    fig, ax = plt.subplots(figsize=(8, 8))

    def init():
        ax.clear()
        ax.set_axis_off()  

    def update(i):
        ax.clear()  
        G = latt_list[i]
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=10, node_color='skyblue', alpha=0.6)
        edge_widths = [d['diameter'] * 3.0 for _, _, d in G.edges(data=True)]
        set1 = set(inf_edges)
        set2 = set(G.edges())
        #set3 = set(merged_edges)
        inf = list(set1-set2)
        #merged = list(set3-set2)
        nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, alpha=1.)
        nx.draw_networkx_edges(G,pos,ax=ax,edge_color='blue',style='dashed',edgelist=inf)
        #nx.draw_networkx_edges(G,pos,ax=ax,edge_color='red',style='dashed',edgelist=merged)
        
        # Draw node labels (if desired)
        # nx.draw_networkx_labels(G, pos, ax=ax, font_size=10)

    # Create the animation, updating the graph every frame
    ani = FuncAnimation(fig, update, frames=len(latt_list), init_func=init, repeat=False, interval=500)

    # Save the animation as a GIF
    ani.save(fp, writer="pillow", fps=1)

def draw(G,pos,inf_edges,i):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.clear()
    ax.set_axis_off()    
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=10, node_color='skyblue', alpha=0.6)
    edge_widths = [d['diameter'] * 3.0 for _, _, d in G.edges(data=True)]
    set1 = set(inf_edges)
    set2 = set(G.edges())
    inf = list(set1-set2)
    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, alpha=1.)
    nx.draw_networkx_edges(G,pos,ax=ax,edge_color='blue',style='dashed',edgelist=inf)
    fig.savefig(f"first_stage_{i}")

edges,nodes,pos,neighbors  = lattice(20,20)


d_init = 0.1
pos_for_export=[]
for p in pos:
    pos_for_export.append(p[1])
'''
save_to_file(edges,"edges.txt")
save_to_file(pos_for_export,"pos.txt")
save_to_file(neighbors,"neighbors.txt")
save_to_file([d_init for i in range(len(edges))],"D.txt")

subprocess.run(["g++", "-o", "computation_loop2.0", "computation_loop2.0.cpp"], check=True)
start_time = time.perf_counter()

subprocess.run(["./computation_loop2.0", "16000", "0.5"], capture_output=False, check=True)

end_time = time.perf_counter()

execution_time = end_time - start_time
print(f"Cpp computation time: {execution_time} seconds")
'''
timeline = import_timeline()
inf_edges = import_edges("infinite_edges.txt")
merged_edges= import_edges("merged_edges.txt")
print(len(timeline))
start_time = time.perf_counter()
for i in {0,73,109,146}:
    draw(timeline[i],pos_for_export,inf_edges,i)

end_time = time.perf_counter()

execution_time = end_time - start_time
print(f"Gif generation time: {execution_time} seconds")
'''
fig, ax = plt.subplots()
G = graphs[-1]
nx.draw_networkx(G,pos=pos_for_export,ax=ax,node_size=5)#uncomment to add visual widths
plt.show()'''