# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 18:42:50 2023

@author: Alex
"""
import os, json, itertools
from ase.io import read
from ase import Atoms
import numpy as np

def readin(fname):
    f = open(fname, 'r')
    content = f.read()
    return content

class pvchem:
    def load_mol(self, filename):
        self.mol = read(filename)
        self.connections = find_connections(self.mol)
        self.mol.positions -= self.mol.positions.min(axis=0)
        
    def make_bond(self, conn):
        print("conn a0:", conn["a0"])
        a0 = self.mol[conn["a0"]]
        a1 = self.mol[conn["a1"]]
        radius = 0.1
        cylinder = ""
        cylinder += "cylinder {\n"
        cylinder += f"<{a0.position[0]}, {a0.position[1]}, {a0.position[2]}>, "
        cylinder += f"<{a1.position[0]}, {a1.position[1]}, {a1.position[2]}>, {radius}\n"
        cylinder += "pigment { rgbt <0.75, 0.75, 0.75, 0> }\n"
        cylinder += "}\n"
        return cylinder
        
    def write(self, fname):
        with open(fname, 'w') as povout:
            povout.write(self.defaults)
            for i in range(len(self.mol)):
                povout.write("\n\n")
                atom = self.mol[i]
                povout.write("sphere {\n")
                povout.write("    <{:.4f}, {:.4f}, {:.4f}>, 0.33\n".format(*atom.position))
                povout.write("    pigment { rgbt <0.75, 0.75, 0.75, 0> }\n")#.format(atom.position[0]))
                povout.write("}\n")
            
            for connection in self.connections:
                if self.connections[connection]['Type'] == "Bond":
                    povout.write(self.make_bond(self.connections[connection]))
                print(connection, self.connections[connection])
        
    def __init__(self):
        self.defaults = readin(os.path.join(os.path.dirname(__file__), "defaults.pov"))

        
def find_connections(mol):
    # Determine the things we are going to vary
    Connections = {}
    for i in range(len(mol)):
        d = mol.get_distances(i, indices=np.arange(0, len(mol)))
        for j in np.where((d > 0.01) & (d < 1.8))[0].tolist():
            key = f"{i}-{j}"
            Connections[key] = {"Type":"Bond", "a0": i, "a1": j,
                          "s0": mol[i].symbol, "sj": mol[j].symbol,
                          "val": float(d[j])}
    for triplet in itertools.permutations(np.arange(len(mol)).tolist(), 3):
        i,j,k = triplet
        if mol.get_distance(i, j) > 1.8 or mol.get_distance(j, k) > 1.8:# or mol.get_distance(j, k) > 3.0:
            continue
        key=f"{i}-{j}-{k}"
        angle = float(mol.get_angle(*triplet))
        Connections[key] = {"Type":"Angle", "a0": i, "a1": j, "a2": k,
                       "s0": mol[i].symbol, "s1": mol[j].symbol, "s2": mol[k].symbol,
                       "val": angle}
    for quadruplet in itertools.permutations(np.arange(len(mol)).tolist(), 4):
        i,j,k,l = quadruplet
        if mol.get_distance(i, j) > 1.8 or mol.get_distance(j, k) > 1.8 or mol.get_distance(k,l) > 1.8:
            continue
        key=f"{i}-{j}-{k}-{l}"
        dihedral = float(mol.get_dihedral(*quadruplet))
        Connections[key] = {"Type":"Dihedral", "a0": i, "a1": j, "a2": k, "a3": l,
                       "s0": mol[i].symbol, "s1": mol[j].symbol, "s2": mol[k].symbol, "s3": mol[l].symbol,
                       "val": dihedral}
    with open(f"Connections.json", 'w') as jout:
        json.dump(Connections, jout, indent=4)
    return Connections