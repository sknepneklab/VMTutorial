
# \file Edge.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 04-Dec-2024
# Hadles edges in the half-edge implementation of the mesh 
#

class Edge:

    def __init__(self,idx):
        self.idx = idx 
        self.vi = None 
        self.vj = None 
        self.he = None