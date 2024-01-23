# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 18:42:50 2023

@author: Alex
"""
import os, json, itertools, pandas
from ase.io import read
from ase import Atoms
import numpy as np
from colour import Color

bond_cutoffs = pandas.DataFrame()
bond_cutoffs.at["H", "H"] = 1.1
bond_cutoffs.at["H", "C"] = 1.4
bond_cutoffs.at["H", "N"] = 1.4
bond_cutoffs.at["H", "O"] = 1.5
bond_cutoffs.at["H", "Cl"] = 1.4
bond_cutoffs.at["C", "C"] = 1.7
bond_cutoffs.at["C", "N"] = 1.7
bond_cutoffs.at["C", "O"] = 1.7
bond_cutoffs.at["C", "Cl"] = 1.9
bond_cutoffs.at["N", "N"] = 1.7
bond_cutoffs.at["N", "O"] = 1.7
bond_cutoffs.at["N", "Cl"] = 1.9
bond_cutoffs.at["O", "O"] = 1.7
bond_cutoffs.at["O", "Cl"] = 1.9
bond_cutoffs.at["Cl", "Cl"] = 2.0

bond_cutoffs.at["H", "K"] = 0.0
bond_cutoffs.at["C", "K"] = 0.0
bond_cutoffs.at["N", "K"] = 0.0
bond_cutoffs.at["O", "K"] = 0.0
bond_cutoffs.at["Cl", "K"] = 0.0

# CG
bond_cutoffs.at["B", "B"] = 4.75
bond_cutoffs.at["B", "S"] = 0.5
bond_cutoffs.at["S", "S"] = 0.5


for i in bond_cutoffs.index:
    for j in bond_cutoffs.columns:
        bond_cutoffs.at[j,i] = bond_cutoffs.at[i,j]
        
def readin(fname):
    f = open(fname, 'r')
    content = f.read()
    return content

atomic_colours = {"H": [0.75, 0.75, 0.75],
                  "C": Color("grey").get_rgb(),
                  "N": Color("blue").get_rgb(),
                  "K": Color("purple").get_rgb(),
                  "B": Color("pink").get_rgb(),
                  "S": Color("yellow").get_rgb(),
                  }
atomic_radii = {"H": 0.33, "C": 0.456, "N": 0.456, "Cl": 0.525, "K": 1.0, "B": 1.0, "S":1.0}

class pvchem:
    def load_mol(self, filename):
        mol = read(filename)
        self._load(mol)

    def _load(self, mol):
        self.mol = mol
        self.connections = find_connections(self.mol)
        self.mol.positions -= self.mol.positions.min(axis=0)
    def make_bond(self, conn):
        # 2 cylinders per bond, each going to the half way point with their own colour
        a0 = self.mol[conn["a0"]]
        a1 = self.mol[conn["a1"]]
        
        halfway = a1.position + ((a0.position-a1.position)/2)
        colour = atomic_colours[conn["s0"]]
        
        radius = 0.1
        cylinder = ""
        cylinder += "cylinder {\n"
        cylinder += f"<{a0.position[0]}, {a0.position[1]}, {a0.position[2]}>, "
        cylinder += f"<{halfway[0]}, {halfway[1]}, {halfway[2]}>, {radius}\n"
        cylinder += "pigment { rgbt <"+", ".join([str(x) for x in colour]) + ", 0> }\n"
        cylinder += "}\n"
        return cylinder
        
    def write(self, fname, colour=None):
        with open(fname, 'w') as povout:
            povout.write(self.defaults)
            povout.write(f"declare {fname.split('.')[0]} = union "+"{")
            for i in range(len(self.mol)):
                povout.write("\n\n")
                atom = self.mol[i]
                povout.write("sphere {\n")
                povout.write("    <{:.4f}, {:.4f}, {:.4f}>, {:.4f}\n".format(*atom.position, atomic_radii[atom.symbol]))
                if colour is None:
                    povout.write("    pigment { rgbt <"+", ".join([str(x) for x in atomic_colours[self.mol[i].symbol]]) +", 0> }\n")#.format(atom.position[0]))
                else:
                    povout.write("    pigment { rgbt <"+", ".join([str(x) for x in Color(colour).get_rgb()]) +", 0> }\n")#.format(atom.position[0]))
                    
                povout.write("}\n")
            
            print(self.connections)
            for connection in self.connections:
                if self.connections[connection]['Type'] == "Bond":
                    povout.write(self.make_bond(self.connections[connection]))
            povout.write("} //end of declare union\n\n")
            povout.write("//object {"+f" {fname.split('.')[0]} "+"}")
        
    def __init__(self):
        self.defaults = readin(os.path.join(os.path.dirname(__file__), "defaults.pov"))

        
def find_connections(mol):
    # Determine the things we are going to vary
    Connections = {}
    for i in range(len(mol)):
        print( mol[i].symbol, bond_cutoffs.max().max())
        d = mol.get_distances(i, indices=np.arange(0, len(mol)))
        for j in np.where((d > 0.01) & (d < bond_cutoffs.max().max()))[0].tolist():
            cutoff = bond_cutoffs.at[mol[i].symbol, mol[j].symbol]
            if d[j] > cutoff:
                continue
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