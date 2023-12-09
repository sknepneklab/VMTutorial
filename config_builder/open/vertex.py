class Vertex:

    def __init__(self, id, r):
        self.id = id
        self.r = r
        self.type = "regular"
        self.boundary = False
        self.neigh = []
        self.faces = []
        self.constraint = "none"
