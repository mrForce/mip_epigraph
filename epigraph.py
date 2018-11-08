import collections
from pyscipopt import Model, quicksum, quickprod

"""
I will begin by taking the proteins, cutting them into k-mers.


This takes in a list of proteins, and a number k, and returns a collections.Counter object, mapping each epitope to the # of times it appears in the proteins
"""
def cut_proteins(proteins, k):
    epitopes = collections.Counter()
    for protein in proteins:
        for i in range(0, len(protein) - k):
            epitope = protein[i:(i + k)]
            epitopes[epitope] += 1
    return epitopes


def solve(epitopes):
    edge_variables = []
    model = Model()
    for epitope,count in epitopes.items():
        
