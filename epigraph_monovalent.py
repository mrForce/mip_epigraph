import collections
import argparse
from Bio import SeqIO
from pyscipopt import Model, quicksum, quickprod



"""
I will begin by taking the proteins, cutting them into k-mers.


This takes in a list of proteins, and a number k, and returns a collections.Counter object, mapping each epitope to the # of times it appears in the proteins
"""
def cut_proteins(protein_location, k):
    epitopes = collections.Counter()
    with open(protein_location, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            protein = str(record.seq)
            for i in range(0, len(protein) - k + 1):
                epitope = protein[i:(i + k)]
                epitopes[epitope] += 1
    return epitopes


def solve(epitopes):
    """
    Map each vertex to the set of incoming edge variables. 
    """
    incoming_edges = collections.defaultdict(list)
    outgoing_edges = collections.defaultdict(list)    
    source_edges = []
    sink_edges = []
    model = Model()
    print('got the model')
    print('size of epitopes: ' + str(len(epitopes)))
    for u, u_count in epitopes.items():
        #add incoming edge from source node
        variable = model.addVar(vtype='B')
        incoming_edges[u].append((variable, u))
        source_edges.append((u, variable))
        #add outgoing edge to sink node
        variable = model.addVar(vtype='B')
        outgoing_edges[u].append((variable, None))
        sink_edges.append(variable)
        for v, v_count in epitopes.items():
            if (not u == v) and u[1::] == v[0:-1]:
                #add a variable for the edge.
                variable = model.addVar(vtype='B', lb=0, ub=1)
                incoming_edges[v].append((variable, u))
                outgoing_edges[u].append((variable, v))
    
    print('added variables')
    
    for u, u_count in epitopes.items():
        #for each vertex, make sure at most one incoming edge is part of the path
        model.addCons(quicksum([x[0] for x in incoming_edges[u]]) <= 1)
        #balance the outgoing and incoming edges for each node
        model.addCons(quicksum([x[0] for x in incoming_edges[u]]) == quicksum([x[0] for x in outgoing_edges[u]]))
    #exactly one edge coming out of the source is part of the path.
    model.addCons(quicksum([x[1] for x in source_edges]) == 1)
    #exactly one edge going into the sink is part of the path
    model.addCons(quicksum(sink_edges) == 1)
    print('added constraints')
    #now, formulate the objective function.
    model.setObjective(quicksum([u_count*quicksum([x[0] for x in incoming_edges[u]]) for u,u_count in epitopes.items()]), 'maximize')
    print('set objective')
    model.optimize()
    starting_vertex = None
    for source_edge in source_edges:
        #sometimes it's slightly off by a tiny bit (like 10^-12 or something)
        if round(model.getVal(source_edge[1])) == 1:
            starting_vertex = source_edge[0]
            break
    assert(starting_vertex)
    path = [starting_vertex]
    while True:
        print(starting_vertex)
        edges = outgoing_edges[starting_vertex]
        variable = None
        next_vertex = None
        for var, seq in edges:
            if round(model.getVal(var)) == 1:
                variable = var
                next_vertex = seq
                break
        print('variable')
        print(variable)
        
        if next_vertex:
            path.append(next_vertex)
            starting_vertex = next_vertex
        else:
            print('breaking')
            break
    print('path')
    print(path)
    assert(len(path) == len(set(path)))
    synthetic = path[0] + ''.join(x[-1] for x in path[1::])
    
    return synthetic
        
                
parser = argparse.ArgumentParser()
parser.add_argument('proteins', help='proteins.fasta file')
parser.add_argument('k', type=int, help='length')
args = parser.parse_args()


epitopes = cut_proteins(args.proteins, args.k)

synthetic = solve(epitopes)
print('synthetic')
print(synthetic)
