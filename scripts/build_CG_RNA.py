import sys
import numpy as np
from HyresBuilder import RNABuilder

pdb = sys.argv[1]
seq = sys.argv[2]
RNABuilder.build(seq, pdb)
