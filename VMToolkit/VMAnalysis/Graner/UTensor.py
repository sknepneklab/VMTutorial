
# \file UTensor.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 09-Dec-2023
# Computes the U tensor (Graner statistics) and produces a VTP file for it
#

import json
import os

import numpy as np

from ..utils.HalfEdge import *
from .Tensor import Tensor
from .MTensor import MTensor


class UTensor:

    def __init__(self, frame, ref):
        ext_1 = os.path.splitext(frame)[1]
        ext_2 = os.path.splitext(ref)[1]
        if (ext_1 != '.json' and ext_1 != '.bz2') or (ext_2 != '.json' and ext_2 != '.bz2'):
            raise Exception(
                'JSON files with two consecutive frames of the mesh state have to be provided.')
        self.M0 = MTensor(ref)
        self.M = MTensor(frame)
        self.U = Tensor(self.M.M.N)
        for i in range(self.M.M.N):
            if (self.M.mesh.num_inner_faces != self.M0.mesh.num_inner_faces) and len(self.M.original_face_id) > 0:
                ref_idx = self.M.original_face_id[i]
            else:
                ref_idx = i
            ev_ref, R_ref = np.linalg.eigh(self.M0.M.T[ref_idx, :, :])
            ev_frame, R_frame = np.linalg.eigh(self.M.M.T[i, :, :])
            log_ref = np.diag(np.log(ev_ref))
            log_frame = np.diag(np.log(ev_frame))
            logM0 = np.dot(np.linalg.inv(R_ref), np.dot(log_ref, R_ref))
            logM = np.dot(np.linalg.inv(R_frame), np.dot(log_frame, R_frame))
            self.U.T[i, :, :] = 0.5*(logM - logM0)

    def plot_vtk_tensor(self, filename):
        self.U.plot_vtk_tensor(filename, self.M.mesh, 'U_tensor')
