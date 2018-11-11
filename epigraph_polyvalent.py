import collections
import argparse
from Bio import SeqIO
from pyscipopt import Model, quicksum, quickprod
from Bio.Alphabet.IUPAC import IUPACProtein
import itertools

"""
I will begin by taking the proteins, cutting them into k-mers.


This takes in a list of proteins, and a number k, and returns a collections.Counter object, mapping each epitope to the # of times it appears in the proteins
"""
def cut_proteins(protein_location, k):
    epitopes = collections.Counter()
    with open(protein_location, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            protein = ''.join(list(filter(lambda x: x in IUPACProtein.letters, str(record.seq))))
            for i in range(0, len(protein) - k + 1):
                epitope = protein[i:(i + k)]
                epitopes[epitope] += 1
    return epitopes


def solve(epitopes, num_synthetic_antigens):
    """
    
    If there are n synthetic antigens, each edge is cloned across the n synthetic antigens.
    
    """
    edge_variables = []
    incoming_edges = collections.defaultdict(lambda: collections.defaultdict(list))
    outgoing_edges = collections.defaultdict(lambda: collections.defaultdict(list))    
    source_edges = collections.defaultdict(list)
    sink_edges = collections.defaultdict(list)
    z_variables = dict()
    model = Model()
    print('got the model')
    print('size of epitopes: ' + str(len(epitopes)))
    print('number of synthetic antigens to generate')
    print(num_synthetic_antigens)
    for u, u_count in epitopes.items():
        z_variables[u] = model.addVar(vtype='B')
    for u, u_count in epitopes.items():
        for k in range(0, num_synthetic_antigens):
            #add incoming edge from source node
            variable = model.addVar(vtype='B')
            incoming_edges[k][u].append((variable, None))
            source_edges[k].append((u, variable))
            #add outgoing edge to sink node
            variable = model.addVar(vtype='B')
            outgoing_edges[k][u].append((variable, None))
            sink_edges[k].append(variable)
            for v, v_count in epitopes.items():
                if (not u == v) and u[1::] == v[0:-1]:
                    #add a variable for the edge.
                    variable = model.addVar(vtype='B', lb=0, ub=1)
                    edge_variables.append(variable)
                    #edge goes from u to v. 
                    model.addCons(variable <= z_variables[v])
                    incoming_edges[k][v].append((variable, u))                    
                    outgoing_edges[k][u].append((variable, v))
    
    
    for v, v_count in epitopes.items():
        #print('list')
        #print(list(itertools.chain.from_iterable([[y[0] for y in x[v]] for x in incoming_edges.values()])))
        model.addCons(quicksum(list(itertools.chain.from_iterable([[y[0] for y in x[v]] for x in incoming_edges.values()]))) >= z_variables[v])

    for k in range(0, num_synthetic_antigens):
        model.addCons(quicksum([x[1] for x in source_edges[k]]) == 1)
        #exactly one edge going into the sink is part of the path
        model.addCons(quicksum(sink_edges[k]) == 1)
        for u, u_count in epitopes.items():
            #for each vertex, make sure at most one incoming edge is part of the path
            model.addCons(quicksum([x[0] for x in incoming_edges[k][u]]) <= 1)
            #balance the outgoing and incoming edges for each node
            model.addCons(quicksum([x[0] for x in incoming_edges[k][u]]) == quicksum([x[0] for x in outgoing_edges[k][u]]))
    print('added constraints')
    #now, formulate the objective function.
    model.setObjective(quicksum([u_count*z_variables[u] for u,u_count in epitopes.items()]), 'maximize')
    print('set objective')
    model.optimize()
    total = 0
    num_true = 0
    for var in edge_variables:
        value = model.getVal(var)
        total += value
        if value >= 0.0000000001:
            num_true += 1
        
    print('total')
    print(total)
    print('num true')
    print(num_true)
    pool = []
    for k in range(0, num_synthetic_antigens):
        starting_vertex = None                
        for source_edge in source_edges[k]:
            #sometimes it's slightlya off by a tiny bit (like 10^-12 or something)
            if round(model.getVal(source_edge[1])) == 1:
                starting_vertex = source_edge[0]
                break
        assert(starting_vertex)
        path = [starting_vertex]
        while True:
            edges = outgoing_edges[k][starting_vertex]
            variable = None
            next_vertex = None
            for var, seq in edges:
                if round(model.getVal(var)) == 1:
                    variable = var
                    next_vertex = seq
                    break
        
            if next_vertex:
                path.append(next_vertex)
                starting_vertex = next_vertex
            else:
                print('breaking')
                break
        assert(len(path) == len(set(path)))
        synthetic = path[0] + ''.join(x[-1] for x in path[1::])
        pool.append(synthetic)
    return pool
        
                
parser = argparse.ArgumentParser()
parser.add_argument('proteins', help='proteins.fasta file')
parser.add_argument('k', type=int, help='length')
parser.add_argument('n', type=int, help='num synthetic antigens')
args = parser.parse_args()


epitopes = cut_proteins(args.proteins, args.k)

pool = solve(epitopes, args.n)
print('synthetic')
for synthetic in pool:
    print(synthetic)
