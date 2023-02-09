#! /Users/liangying/Documents/Code/FreeEnergy/.venv/bin/python
import numpy as np
from lib import vector_displacement

poscar_distort = '/Users/liangying/Documents/Server/T3/Work/TiSe2/bulk-2x2x2/ISIF2-GGA-sym-highk-encut-700/4step/0-Relax/CONTCAR'
poscar_normal = '/Users/liangying/Documents/Server/T3/Work/TiSe2/bulk-1x1x1/13step/1-Sc/POSCAR'

displacement = vector_displacement(poscar_normal, poscar_distort)
