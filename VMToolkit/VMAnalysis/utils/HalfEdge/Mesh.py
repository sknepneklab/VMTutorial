#
# \file Mesh.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 04-Dec-2024
# Handles the mesh in the half-edge implementation of the mesh

from .Box import Box
from .Vec import Vec
from .Vertex import Vertex 
from .Edge import Edge 
from .HalfEdge import HalfEdge
from .Face import Face 
import json
import bz2
import numpy as np
import pathlib
try:
    import vtk
except ImportError:
    raise Exception('Python VTK support needs to be installed. Please check VTK documentation (vtk.org) on how to install it.')


class Mesh:

    def __init__(self):
        self.vertices = []
        self.edges = []
        self.he = []
        self.faces = []
        self.hep_map = {}     # Dictionary of pairs (stack of temporary boundary edges)
        self.box = None
        self.system = {}

    def __add_vertex(self, idx, r):
        self.vertices.append(Vertex(idx,r))
    
    def __add_face(self, idx, verts, erased = False):
        face = Face(idx)    #  Create new face
        face.erased = erased 
        if not erased:
            he_list = []        #  Local storage for prev/next bookkeeping 
            N = len(verts)
            for i in range(N):  # Loop over all vertices in the face
                vi = self.vertices[verts[i]]              #  First vertex
                vj = self.vertices[verts[(i+1) % N]]      #  and its follower, that can be wrapped around
                new_pair = False           #  flag that checks if we created a new pair of half-edges 
                if not (vi,vj) in self.hep_map:    #  if the half-edge does not exist
                    he = HalfEdge(len(self.he))    #  create it
                    hep = HalfEdge(he.idx+1)       #  and set its pair
                    he.pair = hep                  #  set them as pairs of each other
                    hep.pair = he
                    hep.boundary = True            #  per construction, pair is boundary 
                    self.hep_map[(vj,vi)] = hep
                    new_pair = True                #  set new_pair flag to True for post-processing
                    edge = Edge(len(self.edges))
                    edge.vi = vi 
                    edge.vj = vj
                    edge.he = he
                    self.edges.append(edge)
                else:                              #  if the half-edge exists retrieve it
                    he = self.hep_map[(vi,vj)]
                    he.boundary = False            #  per construction, it cannat be boundary 
                    self.hep_map.pop((vi,vj), None)    #  remove it from the stack
                # Do the bookkeeping of the connectivity 
                he.vfrom = vi 
                he.vto = vj 
                he.face = face
                he_list.append(he)
                if i > 0:
                    he.prev = he_list[i-1] 
                    he_list[i-1].next = he
                if i == N-1:
                    he.next = he_list[0]
                    he_list[0].prev = he
                    face.he = he_list[0]  # This makes sure that we can read off things such as myosin and tension
                self.vertices[vi.idx].he = he
                # Add new pair to the list of all half-edges
                if new_pair:
                    self.he.append(he)
                    self.he.append(hep)
        self.faces.append(face)

    def __build_boundary(self):
        if len(self.hep_map) == 0:
            raise RuntimeError('build_boundary() needs to be called after all faces have been added.')
        # After all inner faces have been added, what is left in the hep_map dictionary
        # are boundary edges. 
        # We need to connect them to each other in the proper order
        face = Face(-1)
        face.outer = True
        for (vi,vj) in self.hep_map.keys():
            he = self.hep_map[(vi,vj)]
            he.vfrom = vi 
            he.vto = vj 
            he.face = face
            hev = vj.he
            while not hev.boundary:
                hev = hev.prev.pair
            he.next = hev 
            hev.prev = he
        face.he = he
        self.outer = face
        self.outer.he = he   # This is the fictitious outer face. 
        self.faces.append(face)

    def read(self, f1, f2 = None, box = None, include_erased = True):
        if pathlib.Path(f1).suffix.lower() == '.json':
            self.read_json(f1, include_erased = include_erased)
        else:
            self.read_data(f1, f2, box)

    def read_data(self, f1, f2, box, A0):
        with open(box, 'r') as box_file:
            lx, ly = list(map(float, box_file.read().strip().split()))
            self.box = Box([lx,0], [0,ly])
        with open(f1, 'r') as vert_file:
            for line in vert_file:
                sp_line = line.strip().split()
                vid = int(sp_line[0])
                r = [float(sp_line[1]), float(sp_line[2])]
                self.__add_vertex(vid, Vec(r, self.box))
                self.vertices[-1].type = 'regular'
                self.vertices[-1].boundary = False
                self.vertices[-1].constraint = 'none'
                if len(sp_line) > 3:
                    if int(sp_line[3]) == 1:
                        self.vertices[-1].constraint = 'fixed'
                self.vertices[-1].erased = False
        with open(f2, 'r') as face_file:
            for line in face_file:
                sp_line = line.strip().split()
                fid = int(sp_line[0])
                num_vert = int(sp_line[1])
                verts = list(map(int, sp_line[2:num_vert+2]))
                n = list(map(float, sp_line[num_vert+2:-1]))
                prop = float(sp_line[-1]) 
                self.__add_face(fid, verts) 
                self.faces[-1].A0 = A0
                self.faces[-1].n = n
                self.faces[-1].type = "passive"
                self.faces[-1].property = prop
        self.num_inner_faces = list(map(lambda x: not x.outer, self.faces)).count(True)
    
    def read_json(self,meshfile,include_erased=True):
        with (open(meshfile,'r') if meshfile.endswith(".json") else open(meshfile,'rb') ) as mf:
            if meshfile.endswith("bz2") :
                json_data = bz2.decompress(mf.read())                
            else :
                json_data = mf.read()
            mesh = json.loads(json_data)
            if 'system' in mesh: #loads metadata
                for (key, value) in mesh['system'].items():
                    self.system[key] = value

            if 'box' in mesh["mesh"] and mesh["mesh"]["box"]["periodic"]:
                a, b = [], []
                if 'lx' in mesh["mesh"]["box"]:
                    a = [mesh["mesh"]["box"]["lx"], 0.0]
                    if 'ly' in mesh["mesh"]["box"]:
                        b = [0.0, mesh["mesh"]["box"]["ly"]]
                    else:
                        b = [0.0, mesh["mesh"]["box"]["lx"]]
                if 'a' in mesh["mesh"]["box"] and 'b' in mesh["mesh"]["box"]:
                    a = mesh["mesh"]["box"]["a"]
                    b = mesh["mesh"]["box"]["b"]
                self.box = Box(a,b)
            else:
                self.box = None
            vid = 0
            vidmap = {}
            for v in mesh["mesh"]["vertices"]:
                vidmap[v["id"]] = vid
                if include_erased or not v["erased"]:
                    self.__add_vertex(vid, Vec(v["r"],self.box))
                    self.vertices[-1].type = v["type"]
                    self.vertices[-1].boundary = v["boundary"]
                    self.vertices[-1].constraint = v["constraint"]
                    self.vertices[-1].erased = v["erased"]
                    vid += 1
            has_outer = False
            for f in mesh["mesh"]["faces"]:
                if not f["outer"]:
                    fverts = [vidmap[vid] for vid in f['vertices']]
                    if "erased" in f:
                        self.__add_face(f["id"], fverts, f["erased"])
                    else:
                        self.__add_face(f["id"], fverts)
                    if "type" in f:
                        self.faces[-1].type = f["type"]
                    if ("myo" in f or "tension" in f) and not f["erased"]:           # Read in myosin and tension if they exist
                        face = self.faces[-1]
                        he = face.he 
                        first = face.he 
                        i = 0
                        while True:
                            if "myo" in f:
                                he.myo = f["myo"][i]
                            if "tension" in f:
                                he.tension = f["tension"][i]
                            i += 1
                            he = he.next 
                            if he.idx == first.idx:
                                break
                    if "kappa" in f:
                        self.faces[-1].params["kappa"] = f["kappa"]
                    if "gamma" in f:
                        self.faces[-1].params["gamma"] = f["gamma"]
                    if "lambda" in f:
                        self.faces[-1].params["lambda"] = f["lambda"]
                    if "A0" in f:
                        self.faces[-1].params["A0"] = f["A0"]
                        self.faces[-1].A0 = f["A0"]
                    if "original_id" in f:
                        self.faces[-1].params["original_id"] = f["original_id"]
                else:
                    has_outer = True
            if has_outer:
                self.__build_boundary()
        
            self.num_inner_faces = list(map(lambda x: not x.outer, self.faces)).count(True)
    
    def dat_to_json(self, vert_file, cell_file, box_file, json_file, A0 = 1):
        self.read_data(vert_file, cell_file, box_file, A0)
        with open(json_file, 'w') as out:
            json.dump(self.to_json(), out, indent=4)

    def displace(self, dr):
        if dr.size != 2*len(self.vertices):
            raise Exception('Displacement vector has to have length equal to twice the number of vertices.')
        i = 0
        for v in self.vertices:
            v.r += Vec([dr[i], dr[i+1]])
            i += 2

    def flip_edge(self, eid):
        if eid < 0 or eid >= len(self.edges):
            raise ValueError('Wrong edge id.')
        # Get Edge
        edge = self.edges[eid]
        he = edge.he
        hep = he.pair
        vi = edge.vi
        vj = edge.vj
        # Rotate vertices by pi/2 around edge centre
        rc = Vec(0.5*(vi.r.r + vj.r.r), self.box)
        ric = vi.r - rc
        ric.rotate(0.5*np.pi)
        rjc = vj.r - rc
        rjc.rotate(0.5*np.pi)
        vi.r = ric + rc
        vj.r = rjc + rc
        # Update connectivity 
        tmp_he_1 = he.prev.pair 
        tmp_he_2 = hep.prev.pair

        hep.next.pair.next = he 
        he.next.pair.next = hep 

        hep.prev.pair.prev = hep 
        he.prev.pair.prev = he 

        he.prev.next = he.next 
        he.next.prev = he.prev 

        hep.prev.next = hep.next 
        hep.next.prev = hep.prev

        he.prev = hep.next.pair 
        hep.prev = he.next.pair 
        he.next = tmp_he_1
        hep.next = tmp_he_2

        # update vertices
        he.next.vfrom = vj 
        hep.next.vfrom = vi

        # update face
        he.face = he.prev.face
        hep.face = hep.prev.face

    def nneighbours(self, idx, numneigh = 2):
        nn = numneigh
        nb = [idx]
        while nn > 0:
          shell = []
          for i in nb:
            for n in self.faces[i].neighbours():
              if not n in shell:
                shell.append(n)
          for s in shell:
            if not s in nb:
              nb.append(s)
          nn -= 1
        return nb
    
    def get_original_face(self, id):
        for f in self.faces:
            if 'original_id' not in f.params:
                raise Exception('get_original_face can only be used in original_if flag has been set.')
            if f.params['original_id'] == id:
                return f 
        else:
            return None
    
    def scale(self, scale, A = 0):
        a = scale*self.box.a
        b = scale*self.box.b 
        self.box = Box(a.tolist(), b.tolist())
        for v in self.vertices:
            dr = np.random.normal(scale=A,size=2)
            v.r = Vec(scale*v.r.r, self.box)
            v.r = v.r + Vec(dr, self.box)
        for f in self.faces:
            f.A0 = f.area()
            if 'A0' in f.params:
                f.params['A0'] = f.A0
        
    def to_json(self, A0 = None):
        jsonData = {}
        jsonData["mesh"] = {}
        jsonData["mesh"]["vertices"] = []
        if self.box != None:
            jsonData["mesh"]["box"] = {"periodic": True, "lx": self.box.a[0], "ly": self.box.b[1]}
            #jsonData["mesh"]["box"] = {"periodic": True, "a": self.box.a.tolist(), "b": self.box.b.tolist()}
        for v in self.vertices:
            vd = {}
            vd["id"] = v.idx 
            vd["r"] = v.r.to_list()
            vd["type"] = v.type
            vd["erased"] = False 
            vd["boundary"] = v.boundary 
            vd["constraint"] = v.constraint
            jsonData["mesh"]["vertices"].append(vd)
        jsonData["mesh"]["faces"] = []
        for f in self.faces:
            face = {}
            face["id"] = f.idx 
            face["outer"] = f.outer
            face["type"] = f.type
            if f.property:
                face["property"] = f.property 
            verts = []
            he = f.he
            first = f.he
            while True:
                verts.append(he.vfrom.idx)
                he = he.next
                if he.idx == first.idx:
                    break
            face["vertices"] = verts
            for (key, value) in f.params.items():
                face[key] = value
            if A0:
                face["A0"] = A0
            else:
                face["A0"] = f.A0
            jsonData["mesh"]["faces"].append(face)
        return jsonData

    def write_json(self, json_file, A0 = None):
        with open(json_file, 'w') as out:
            json.dump(self.to_json(A0 = A0), out, indent=4)

    def write_vtp(self, fname):
        points = vtk.vtkPoints()
        faces = vtk.vtkCellArray()
        for (i, v) in enumerate(self.vertices):
            X, Y = v.r.r
            points.InsertNextPoint([X, Y, 0.0])

        areas = vtk.vtkDoubleArray()
        areas.SetNumberOfComponents(1)
        areas.SetName("Area")
        perims = vtk.vtkDoubleArray()
        perims.SetNumberOfComponents(1)
        perims.SetName("Perimeter")
        if self.faces[0].property:
            prop = vtk.vtkDoubleArray()
            prop.SetNumberOfComponents(1)
            prop.SetName("Property")
        for f in self.faces:
            add_face = not (f.crosses_boundary() or f.outer)
            if add_face:
                face = vtk.vtkTriangle()
                face.GetPointIds().SetNumberOfIds(f.num_sides())
                he = f.he
                first = f.he
                i = 0
                while True:
                    face.GetPointIds().SetId(i, he.vfrom.idx)
                    he = he.next
                    if he.idx == first.idx:
                        break
                    i += 1
                faces.InsertNextCell(face)
                areas.InsertNextValue(f.area())
                perims.InsertNextValue(f.perimeter())
                if f.property:
                    prop.InsertNextValue(f.property)

        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetPolys(faces)
        polyData.GetCellData().AddArray(areas)
        polyData.GetCellData().AddArray(perims)
        if self.faces[0].property:
            polyData.GetCellData().AddArray(prop)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fname)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(polyData)
        else:
            writer.SetInputData(polyData)
        writer.SetDataModeToAscii()
        writer.Write()    

    def write_vtp_nematic(self, fname, scale = 1.0):
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        i = 0
        for f in self.faces:
            if not f.outer:
                r = f.centroid()
                points.InsertNextPoint([r.r[0], r.r[1], 0.0])
                xi = r.r[0] + 0.5*scale*f.n[0]
                yi = r.r[1] + 0.5*scale*f.n[1]
                points.InsertNextPoint([xi, yi, 0.0])
                l = vtk.vtkLine()
                l.GetPointIds().SetId(0, i)
                l.GetPointIds().SetId(1, i+1)
                lines.InsertNextCell(l)
                xi = r.r[0] - 0.5*scale*f.n[0]
                yi = r.r[1] - 0.5*scale*f.n[1]
                points.InsertNextPoint([xi, yi, 0.0])
                l = vtk.vtkLine()
                l.GetPointIds().SetId(0, i)
                l.GetPointIds().SetId(1, i+2)
                lines.InsertNextCell(l)
                i += 3

        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetLines(lines)
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fname)
        writer.SetInputData(polyData)
        writer.SetDataModeToAscii()
        writer.Write()
