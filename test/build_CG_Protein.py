import sys
import numpy as np
from iConBuilder import ProteinBuilder

pdb = sys.argv[1]
seq = sys.argv[2]
ProteinBuilder.build(seq, pdb)
