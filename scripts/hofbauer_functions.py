import numpy as np
from numpy.linalg import eig
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs


def f(a, x):
    '''
    Our 'tent map' for the dynamical system.

    Parameters
    ----------
    a : float
        The value of the variable in the tent map.
        Can be values between 1 and 2 only.
    x : float
        The current point in the domain of the tent map.
    
    Returns
    -------
    The value of the tent map for specified a and x.
    '''
    return round(a * np.min([x, 1-x]), 5)



def f_L(a, x):
    '''
    The left branch of the tent map.

    Parameters
    ----------
    a : float
        The value of the variable in the tent map.
        Can be values between 1 and 2 only.
    x : float
        The current point in the domain of the tent map.

    Returns
    -------
    The value of the left branch of the tent map for specified a and x.
    '''
    c_0 = 0.5

    if x >= 0 and x <= c_0:
        return round(a * x, 8)
    else:
        return [None]


def f_R(a, x):
    '''
    The right branch of the tent map.

    Parameters
    ----------
    a : float
        The value of the variable in the tent map.
        Can be values between 1 and 2 only.
    x : float
        The current point in the domain of the tent map.

    Returns
    -------
    The value of the left branch of the tent map for specified a and x.
    '''
    c_0 = 0.5

    if x >= c_0 and x <= 1:
        return round(a * (1-x), 8)
    else:
        return [None]
    

def range_f_L(a, interval):
    '''
    Finds the range of values produced by the left branch of the tent map on a given interval.
    NOTE: Finds the range by checking the start and end values of inputted interval.

    Parameters
    ----------
    a : float
        The value of the variable in the tent map.
        Can be values between 1 and 2 only.
    interval : list
        The domain across which to find the range of the branch. 
    '''
    c_0 = 0.5

    if interval[0] > c_0:
        return [None]
    if interval[1] > c_0:
        return [f_L(a, interval[0]), f_L(a, c_0)]
    else:
        return [f_L(a, interval[0]), f_L(a, interval[1])]
    

def range_f_R(a, interval):
    '''
    Finds the range of values produced by the right branch of the tent map on a given interval.
    NOTE: Finds the range by checking the start and end values of inputted interval.

    Parameters
    ----------
    a : float
        The value of the variable in the tent map.
        Can be values between 1 and 2 only.
    interval : list
        The domain across which to find the range of the branch. 
    '''
    c_0 = 0.5

    if interval[1] < c_0:
        return [None]
    if interval[0] < c_0:
        return sorted([f_R(a, c_0), f_R(a, interval[1])])
    else:
        return sorted([f_R(a, interval[0]), f_R(a, interval[1])])
    



def create_digraph(a, limit):
    '''
    Creates and stores a directed graph holding information about Hofbauer extension.

    Parameters
    ----------
    a : float
        The value of the variable in the tent map.
        Can be values between 1 and 2 only.
    limit: int
        The recursion limit when creating the directed graph.
        Corresponds to the number of intervals visited.

    Returns
    -------
    G : DiGraph
        networkx object that stores the infomration about the Hofbauer extension as nodes and labelled, directed edges
    '''
    G = nx.DiGraph()
    recursive_create_node(G, 0, limit, [0.0,1.0], a)
    return G
    

def recursive_create_node(G, counter, limit, interval, a):

    if counter >= limit:
        return # End recursion
    
    interval_tuple = tuple(interval)

    r_l = range_f_L(a, interval) # Compute range of left branch on given interval
    r_l_tuple = tuple(r_l) # Convert list to a tuple for comparison to nodes

    r_r = range_f_R(a, interval) # Compute range of right branch on given interval
    r_r_tuple = tuple(r_r) # Convert list to a tuple for comparison to nodes
 
    if r_l == r_r: # Branches map to the same range
        if r_l_tuple not in G.nodes:
            G.add_node(r_l_tuple) # Create node for this range
        if (interval_tuple, r_l_tuple) not in G.edges:
            G.add_edge(interval_tuple, r_l_tuple, label = 'LR') # Create edge with label indicating both branches map here
        recursive_create_node(G, counter+1, limit, r_l, a) # Recursive call with this range as new interval
    
    else:
        if r_l != [None]: # Within domain of left branch
            if r_l_tuple not in G.nodes: 
                G.add_node(r_l_tuple) # Create node for the range of left branch for this domain
                recursive_create_node(G, counter+1, limit, r_l, a) # Recursive call with this range as new interval
            if (interval_tuple, r_l_tuple) not in G.edges:
                G.add_edge(interval_tuple, r_l_tuple, label = 'L') # Create edge with left label
            
        if r_r != [None]: # Within domain of right branch
            if r_r_tuple not in G.nodes:
                G.add_node(r_r_tuple) # Create node for the range of right branch for this domain
                recursive_create_node(G, counter+1, limit, r_r, a) # Recursive call with this range as new interval
            if (interval_tuple, r_r_tuple) not in G.edges:
                G.add_edge(interval_tuple, r_r_tuple, label = 'R')  # Create edge with right label

def create_adjacency_matricies(a, limit):
    '''
    Takes in a value for a in the tent map function.
    Creates the left, right and combined adjacency matricies for the corresponding directed graph.

    Parameters
    ----------
    a : float
        The value of the variable in the tent map.
        Can be values between 1 and 2 only.
    limit : int
        The recursion limit when creating the directed graph.
        Corresponds to the number of intervals visited.

    Returns
    -------
    nodes_list : list of tuples
        The nodes of the graph, the labels of the columns/rows in the adjacency matricies
    df : np.array
        The labelled adjacency matrix of the DiGraph including edges for both branches.
    df_left : np.array
        The labelled adjacency matrix of the Digraph including edges for the left branch only.
    df_right : np.array
        The labelled adjacncy matrix of the Digraph including edges for the right branch only.
    '''
    G = create_digraph(a, limit)

    nodes_list = list(G.nodes)
    if len(nodes_list) >= 2:
        nodes_list[0], nodes_list[1] = nodes_list[1], nodes_list[0]

    edges_list = []
    for u, v, attr in G.edges(data=True):
        edges_list.append([u,v,attr])

    df = pd.DataFrame(index=nodes_list, columns=nodes_list)
    df.fillna(0, inplace=True) # Create data frame of all zeroes with dimension equal to number of nodes

    for edge in edges_list:
        node1 = edge[0]
        node2 = edge[1]
        label = edge[2]
        if label == {'label': 'LR'}:
            df.at[node2, node1] = 2 # If both branches map we have 2 edges between these nodes
        else:
            df.at[node2, node1] = 1 # Store 1 indicating there is an edge between these nodes

    df_left = pd.DataFrame(index=nodes_list, columns=nodes_list)
    df_right = pd.DataFrame(index=nodes_list, columns=nodes_list)
    df_left.fillna(0, inplace=True)
    df_right.fillna(0, inplace=True)


    for edge in edges_list:
        node1 = edge[0]
        node2 = edge[1]
        label = edge[2]
        if label == {'label': 'L'}:
            df_left.at[node2, node1] = 1
        elif label == {'label': 'R'}:
            df_right.at[node2, node1] = 1
        else:
            df_left.at[node2, node1] = 1
            df_right.at[node2, node1] = 1

    

    return [nodes_list, df.values, df_left.values, df_right.values]


def plot_digraph(a, limit):
    G = create_digraph(a, limit)

    edge_colors = [
        "red" if attr.get("label") == "R" else 
        "green" if attr.get("label") == "L" else 
        "blue" for u, v, attr in G.edges(data=True)
    ]

    pos = nx.shell_layout(G)  # Layout for the graph
    nx.draw(G, pos, with_labels=True,
            node_color="lightblue", edge_color=edge_colors, 
            node_size=2000, font_size=10, font_weight="bold")
    edge_labels = nx.get_edge_attributes(G, "label")  # Extract edge labels

    plt.title(f"DiGraph with a = {a}")
    plt.show()


### Example execution
if __name__ == '__main__':
    matricies = create_adjacency_matricies(1.8, 20)
    print(matricies[1])

plot_digraph(1.8, 50)