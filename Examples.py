# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 18:40:57 2023

@author: Alex
"""

import povray_chem

pv = povray_chem.pvchem()
pv.load_mol("H2.xyz")

pv.write("Test.pov")

#pvengine64 +IH2.pov +OH2.png +W1500 +H1500 +V +D +FN +Q9 +P +UD +UL +UV +A +AM2 +UA

