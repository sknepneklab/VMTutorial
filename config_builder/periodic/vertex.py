class Vertex:

    def __init__(self, id, r):
        self.id = id
        self.r = r
        self.type = "regular"
        self.boundary = False
        self.neigh = []
        self.faces = []
        self.face_centres = []
        self.constraint = "none"
        self.original_idx = id
