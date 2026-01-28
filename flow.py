import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.animation import FuncAnimation
import scipy.sparse.linalg as sla
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse import bmat
from scipy.linalg import qr
import subprocess
import time
import ast
from numpy import where
from scipy.sparse import identity



def position_import(filename):
    pos ={}
    with open(filename, 'r') as f:
        lines = f.readlines()
        for idx, line in enumerate(lines):
            x_str, y_str = line.strip().split()
            x, y = float(x_str), float(y_str)
            pos[idx] = (x, y)
    return pos

def inf_edges_gen(G,all_inf_edges):
    inf_edges = list(set(all_inf_edges)-set(G.edges))
    for e in inf_edges:
        G.add_edge(e[0],e[1],diameter=1e10)
    return sorted(inf_edges)

def import_timeline(path="timeline.txt"):
    graphs = []
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

def length(n1, n2):
    return np.sqrt((n1[0] - n2[0]) ** 2 + (n1[1] - n2[1]) ** 2)


def compute_edge_flow(G,mu,pos,k):
    for u, v in G.edges:
        if G.edges[u, v]['diameter'] == 1e10:
            G.edges[u, v]['flow'] = 'inf'
            continue
        p_u = G.nodes[u]['pressure']
        p_v = G.nodes[v]['pressure']
        d = G.edges[u, v]['diameter']
        l = length(pos[u],pos[v])
        K = (d ** 4) / (128 * mu)
        q = K / l * (p_u - p_v)
        G.edges[u, v]['flow'] = round(k*q,1)
    return G

def fix_diagonal(A):
    n = A.shape[0]
    for i in range(n):
        if A[i, i] == 0:
            A.rows[i] = [i]      
            A.data[i] = [1.0] 
    return A

def inf_man(A,node,b,val,flrow):
    flrow.add(node)  
    A.rows[node] = [node]
    A.data[node] = [1]
    b[node] = val
    return A,b,flrow

def solve_eq(A,b):
    A_dense = A.toarray()
    Q, R = np.linalg.qr(A_dense)
    Q_prim = np.matmul(Q.T,b)
    sol  = np.matmul(np.linalg.inv(R),Q_prim)
    return sol

def gen_C(inf_edges,m,n,frow,lrow,A,b):
    C = lil_matrix((m, n))
    for row_idx, (node0, node1) in enumerate(inf_edges):
    #checking if infinite edge is connected to first row
        if (node0 in frow and node1 in frow) or (node0 in lrow and node0 in lrow):
                _=0
        elif node0 in frow:
            A,b,frow = inf_man(A,node1,b,1,frow)
        elif node1 in frow:
            A,b,frow = inf_man(A,node0,b,1,frow)
        elif node0 in lrow:
            A,b,lrow = inf_man(A,node1,b,0,lrow)
        elif node1 in lrow:
            A,b,lrow = inf_man(A,node0,b,0,lrow)
        C[row_idx, node0] = -1
        C[row_idx, node1] = 1
        if frow & lrow:
            print("The first and the last row are connected via non resistant edges")
            return 0,frow,lrow,A,b
    C_dense = C.toarray()
    rank = np.linalg.matrix_rank(C_dense)
    if rank < C.shape[0]:
        _, _, P = qr(C_dense.T, pivoting=True)  # QR on transpose to find independent rows
        ind_rows = sorted(P[:rank])
        C_cleaned = C[ind_rows, :]
        return C_cleaned,frow,lrow,A,b
    return C,frow,lrow,A,b

def assign_pressures(G,pos,all_nodes,all_inf_edges,mu,frow,lrow):
    
    node_list = sorted(list(G.nodes))
    n = len(all_nodes)
    A = lil_matrix((n, n))
    b = np.zeros(n)
    for node in all_nodes:
        if node in frow:
            A[node, node] = 1
            b[node] = 1
        elif node in lrow: 
            A[node, node] = 1
        elif node in node_list:
            A[node, node] = 0
            for neighbor in G.neighbors(node):
                edge = G.edges[node, neighbor]
                L = length(pos[node],pos[neighbor])
                K = edge.get('diameter') ** 4 / (128 *mu)
                A[node,neighbor] -= K / L
                A[node, node] += K / L
    A = fix_diagonal(A)
    #infinite flow edges
    inf_edges = inf_edges_gen(G,all_inf_edges)
    m = len(inf_edges)
    if m>0:
        C,frow,lrow,A,b = gen_C(inf_edges,m,n,frow,lrow,A,b)
        if type(C)!= lil_matrix :
            return 0
        #creating matrix with Lagrange multipliers
        A_new = bmat([[A, C.T],
              [C, None]], format='csr')
        b_new = np.concatenate([b,np.zeros(A_new.shape[0]-b.shape[0])])
    else:
        A_new = A.tocsr()
        b_new = b
    rank = np.linalg.matrix_rank(A_new.toarray())
    if rank < A_new.shape[0]:
        print(f"WARNING: Matrix A is rank-deficient. Rank={rank}, Shape={A_new.shape}")
    # Solve system
    p = solve_eq(A_new,b_new)
    # Assign pressures back to graph
    for node in G.nodes:
        G.nodes[node]['pressure'] = float(p[node])
    return G

def draw_graph_with_pressures(G, pos,ax,cmap='viridis'):
    pressures = nx.get_node_attributes(G, 'pressure')
    node_colors = [pressures[n] for n in G.nodes]
    nodes = nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=20,
        cmap=plt.get_cmap(cmap),
        ax=ax
    )
    e_label = {(i,j):d['flow'] for i, j, d in G.edges(data=True)}
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color='gray')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=e_label,font_size=6)
    
    #Draw pressure labels
    #pressure_labels = {n: f"{pressures[n]:.2f}" for n in G.nodes}
    #labels={i: j for i,j in enumerate(G.nodes)}
    #nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)

    ax.set_title("Graph with Node Pressures")
    ax.set_aspect('equal')
    return ax

def import_edges(path="infinite_edges.txt"):
    deleted_edges=[]
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in  lines:
            edge_data = ast.literal_eval(line)
            deleted_edges.append((edge_data[0],edge_data[1]))
    return deleted_edges

def flrow_gen(all_nodes,pos):
    frow,lrow = [],[]
    y_values = [pos[node][1] for node in all_nodes]
    min_y = min(y_values)
    max_y = max(y_values)
    for node in all_nodes:
        y = pos[node][1]
        if np.isclose(y, min_y, atol=0.1):
            frow.append(node)   
    for node in all_nodes[-len(frow):]:
        y = pos[node][1]
        if np.isclose(y, max_y, atol= 0.1): 
            lrow.append(node)
    return set(frow),set(lrow)

def gif_gen(latt_list, pos,all_nodes,all_inf_edges,frow,lrow,fp='flogif.gif'):
    fig, ax = plt.subplots(figsize=(15, 15))
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1))  # Set appropriate vmin and vmax
    sm.set_array([])  # Required for colorbar
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("Pressure")
    def init():
        ax.clear()
        ax.set_axis_off()  
    def update(i):
        ax.clear()
        ax.set_axis_off()
        G = latt_list[i]
        p = assign_pressures(G,pos,all_nodes,all_inf_edges,0.89,frow,lrow)
        if p!=0:
            G = p
        else:
            plt.title("INPUT & OUTPUT IS CONNECTED VIA\n NON-RESISTANT EDGES")
        G = compute_edge_flow(G,0.89,pos)
        draw_graph_with_pressures(G, pos,ax,cmap='viridis')
    ani = FuncAnimation(fig, update, frames=len(latt_list), init_func=init, repeat=False, interval=500)
    ani.save(fp, writer="pillow", fps=1)

def draw(G,pos,all_nodes,all_inf_edges,frow,lrow,i):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.clear()
    ax.set_axis_off()
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1))  # Set appropriate vmin and vmax
    sm.set_array([])  # Required for colorbar
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("Pressure")    
    p = assign_pressures(G,pos,all_nodes,all_inf_edges,0.89,frow,lrow)
    if p!=0:
        G = p
    else:
        plt.title("INPUT & OUTPUT IS CONNECTED VIA\n NON-RESISTANT EDGES")
        return 0 
    G = compute_edge_flow(G,0.89,pos,1e6)
    draw_graph_with_pressures(G, pos,ax,cmap='viridis')
    fig.savefig(f"second_stage_{i}")



graphs  = import_timeline()
all_inf_edges = sorted(import_edges("infinite_edges.txt"))
pos = position_import('pos.txt')

print(len(graphs))
all_nodes = sorted(list(graphs[0].nodes))
frow,lrow = flrow_gen(all_nodes,pos)

i=0
for i in {0,int(len(graphs)/2),int(len(graphs)/3),int(len(graphs)/1.5)}: 
    draw(graphs[i],pos,all_nodes,all_inf_edges,frow,lrow,i)

#gif_gen(graphs[-20::5],pos,all_nodes,all_inf_edges,frow,lrow,'flow20x20.gif')