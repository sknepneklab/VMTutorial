try:
    import pyvista as pv
except ImportError:
    raise ImportError('PyVista needs to be installed before tissue_plot can be used.')

import numpy as np 
import matplotlib
from VMToolkit.VM import *

def make_plotter(tissue, cmap = 'hot', colourby = 'area'):
    vertices = []
    faces = []

    box = tissue.vertices()[0].r.box()
    Lx = box[0]
    Ly = box[-1]

    for v in tissue.vertices():
        vert = [v.r.x, v.r.y, 0]
        vertices.append(vert)

    quantity = []
    for (i,c) in enumerate(tissue.cells()):
        verts = [c.neighbours]
        rc = tissue.get_cell_centroid(i)
        omit = False
        for he in FaceCirculator(c):
            vto = he.vto()
            if abs(vto.r.x - rc.x) >= 0.5*Lx or abs(vto.r.y - rc.y) >= 0.5*Ly:
                omit = True 
                break
            verts.append(he.vto().id)
        if not omit:
            faces.append(verts)
            if colourby == 'area':
                quantity.append(tissue.cell_area(i))
            elif colourby == 'perim':
                quantity.append(tissue.cell_perim(i))
            elif colourby == 'type':
                quantity.append(c.type())
            else:
                raise Exception(f'Unknown quantity {colourby}.')
    
    faces = np.hstack(faces)
    quantity = np.array(quantity)

    cm = matplotlib.cm.get_cmap(cmap)


    mesh = pv.PolyData(vertices, faces)
    mesh.cell_data[colourby] = quantity
    
    plotter = pv.Plotter()
    _ = plotter.add_mesh(
        mesh,
        scalars=quantity,
        lighting=False,
        preference='cell',
        show_edges=True,
        line_width = 1,
        scalar_bar_args={'title': colourby},
        show_scalar_bar=True,
        cmap=cm,
    )
    plotter.camera_position = 'xy'
    return plotter



    