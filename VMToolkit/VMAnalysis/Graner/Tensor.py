#
# \file Tensor.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 09-Dec-2023
# \brief Stores and processes a sequence of 2x2 matrices
#

import numpy as np
import vtk


class Tensor:

    def __init__(self, N):
        self.N = N
        self.T = np.zeros((N, 2, 2))
        self.eigvals = np.zeros((N, 2))
        self.eigvecs = np.zeros((N, 2, 2))
        self.type = np.zeros(N, dtype=np.int)
        self.has_eigenvals = False

    def compute_eigvals(self):
        for i in range(self.T.shape[0]):
            e, v = np.linalg.eigh(self.T[i, :, :])
            self.eigvals[i, :] = e
            self.eigvecs[i, :, :] = v
            if np.all(e >= 0):
                self.type[i] = 1
            elif np.all(e < 0):
                self.type[i] = 3
            else:
                self.type[i] = 2
        self.has_eigvals = True

    def mean(self):
        Tavg = np.mean(self.T, axis=0)
        eigval, eigvec = np.linalg.eigh(Tavg)
        vx, vy = eigvec[:, 1]
        alpha = np.arctan2(vy, vx)
        return (eigval, alpha)

    def histogram(self, nbin=25):
        if not self.has_eigenvals:
            self.compute_eigvals()
        anis = self.eigvals[:, 1]/self.eigvals[:, 0]
        disbin = np.linspace(np.min(anis), np.max(anis), nbin + 1)
        anishist, edges = np.histogram(anis, disbin, density=True)
        return anishist, disbin

    def rose(self, nbin=90):
        if not self.has_eigenvals:
            self.compute_eigvals()
        rosebin = np.linspace(-0.5*np.pi, 0.5*np.pi, nbin + 1)
        vx = self.eigvecs[:, 0, 1]
        vy = self.eigvecs[:, 1, 1]
        alpha = np.arctan2(vy, vx)
        hmm = np.where(alpha < -0.5*np.pi)
        alpha[hmm] += np.pi
        hmm2 = np.where(alpha > 0.5*np.pi)
        alpha[hmm2] -= np.pi
        alphahist, edges = np.histogram(alpha, rosebin, normed=True)
        return alphahist, rosebin

    def plot_vtk_tensor(self, filename, mesh, tensorname):
        if not self.has_eigenvals:
            self.compute_eigvals()
        celltypedict = {}
        ct = 1
        for f in mesh.faces:
            if not f.outer:
                if not f.type in celltypedict:
                    celltypedict[f.type] = ct
                    ct += 1
        points = vtk.vtkPoints()
        types = vtk.vtkIntArray()
        celltypes = vtk.vtkIntArray()
        types.SetNumberOfComponents(1)
        types.SetName("EigType")
        celltypes.SetNumberOfComponents(1)
        celltypes.SetName("CellType")
        tensor = vtk.vtkDoubleArray()
        tensor.SetNumberOfComponents(9)
        tensor.SetName(tensorname)
        for f in mesh.faces:
            if not f.outer:
                rc = f.rc(mesh.box)
                points.InsertNextPoint([rc.r[0], rc.r[1], 0.0])
                types.InsertNextValue(self.type[f.idx])
                celltypes.InsertNextValue(celltypedict[f.type])
                tensor.InsertNextTuple([self.T[f.idx, 0, 0], self.T[f.idx, 0, 1],
                                       0.0, self.T[f.idx, 1, 0], self.T[f.idx, 1, 1], 0.0, 0.0, 0.0, 0.0])

        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.GetPointData().AddArray(types)
        polyData.GetPointData().AddArray(celltypes)
        polyData.GetPointData().AddArray(tensor)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(polyData)
        else:
            writer.SetInputData(polyData)
        # writer.SetDataModeToAscii()
        writer.Write()

    def plot_vtk_ellipse(self, filename, mesh, N = 20, scale = 1.0):
        if not self.has_eigenvals:
            self.compute_eigvals()
        
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        t = np.linspace(0, 2*np.pi, N)
        for (idx, f) in enumerate(mesh.faces):
            if not (f.outer or f.erased):
                r = f.rc(mesh.box)     
                ev_1, ev_2 = self.eigvals[idx, 0], self.eigvals[idx, 1]
                evec_1, evec_2 = self.eigvecs[idx, :, 0], self.eigvecs[idx, :, 1] 
                xi = r.r[0] + scale*ev_1*np.cos(t)*evec_1[0] + scale*ev_2*np.sin(t)*evec_2[0]
                yi = r.r[1] + scale*ev_1*np.cos(t)*evec_1[1] + scale*ev_2*np.sin(t)*evec_2[1]
                

                for i in range(N):
                    points.InsertNextPoint([xi[i], yi[i], 0.0])
                    l = vtk.vtkLine()
                    if i > 0:
                        l.GetPointIds().SetId(0, idx*N + i - 1)
                        l.GetPointIds().SetId(1, idx*N + i)
                        lines.InsertNextCell(l)
            
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetLines(lines)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(polyData)
        writer.SetDataModeToAscii()
        writer.Write()

    def plot_vtk_lines(self, filename, mesh, scale = 1.0):
        if not self.has_eigenvals:
            self.compute_eigvals()
        
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        types = vtk.vtkIntArray()
        types.SetNumberOfComponents(1)
        types.SetName("EigSign")
        i = 0
        for (idx, f) in enumerate(mesh.faces):
            if not (f.outer or f.erased):
                r = f.rc(mesh.box)     
                points.InsertNextPoint([r.r[0], r.r[1], 0.0])
                evs = self.eigvals[idx, :]
                evidx = np.argsort(np.abs(evs))       
                ev_1, ev_2 = evs[evidx[0]], evs[evidx[1]]
                evec_1, evec_2 = self.eigvecs[idx, :, evidx[0]], self.eigvecs[idx, :, evidx[1]] 
                xi = r.r[0] + 0.5*scale*np.abs(ev_1)*evec_1[0] 
                yi = r.r[1] + 0.5*scale*np.abs(ev_1)*evec_1[1]
                zi = 0.0
                points.InsertNextPoint([xi, yi, 0.0])
                l = vtk.vtkLine()
                l.GetPointIds().SetId(0, i)
                l.GetPointIds().SetId(1, i+1)
                lines.InsertNextCell(l)
                if ev_1 >=0:
                    types.InsertNextValue(1)
                else:
                    types.InsertNextValue(-1)
                xi = r.r[0] - 0.5*scale*np.abs(ev_1)*evec_1[0] 
                yi = r.r[1] - 0.5*scale*np.abs(ev_1)*evec_1[1]
                zi = 0.0
                points.InsertNextPoint([xi, yi, zi])
                l = vtk.vtkLine()
                l.GetPointIds().SetId(0, i)
                l.GetPointIds().SetId(1, i+2)
                lines.InsertNextCell(l)
                if ev_1 >=0:
                    types.InsertNextValue(1)
                else:
                    types.InsertNextValue(-1)
                xi = r.r[0] + 0.5*scale*np.abs(ev_2)*evec_2[0] 
                yi = r.r[1] + 0.5*scale*np.abs(ev_2)*evec_2[1]
                zi = 0.0
                points.InsertNextPoint([xi, yi, zi])
                l = vtk.vtkLine()
                l.GetPointIds().SetId(0, i)
                l.GetPointIds().SetId(1, i+3)
                lines.InsertNextCell(l)
                if ev_2 >=0:
                    types.InsertNextValue(1)
                else:
                    types.InsertNextValue(-1)
                xi = r.r[0] - 0.5*scale*np.abs(ev_2)*evec_2[0] 
                yi = r.r[1] - 0.5*scale*np.abs(ev_2)*evec_2[1]
                zi = 0.0
                points.InsertNextPoint([xi, yi, zi])
                l = vtk.vtkLine()
                l.GetPointIds().SetId(0, i)
                l.GetPointIds().SetId(1, i+4)
                lines.InsertNextCell(l)
                if ev_2 >=0:
                    types.InsertNextValue(1)
                else:
                    types.InsertNextValue(-1)
                i += 5
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetLines(lines)
        polyData.GetCellData().AddArray(types)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(polyData)
        writer.SetDataModeToAscii()
        writer.Write()
