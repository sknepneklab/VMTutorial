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
    colours = []

    for v in tissue.vertices():
        vert = [v.r.x, v.r.y, 0]
        vertices.append(vert)

    quantity = []
    for (i,c) in enumerate(tissue.cells()):
        if not c.outer:
            verts = [c.neighbours]
            for he in FaceCirculator(c):
                verts.append(he.vto().id)
            faces.append(verts)
            if colourby == 'area':
                quantity.append(tissue.cell_area(i))
            elif colourby == 'perim':
                quantity.append(tissue.cell_perim(i))
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
        cmap=cm
    )
    plotter.camera_position = 'xy'
    return plotter



    