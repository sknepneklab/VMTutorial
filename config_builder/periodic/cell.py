class Cell:

    def __init__(self, id, rc):
        self.id = id
        self.voroid = None
        self.rc = rc
        self.verts = []
        self.type = "passive"
        self.area = None
