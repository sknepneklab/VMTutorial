#
# \file Tensor.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 04-Dec-2024
# Handles operations on vectors in 2d


import numpy as np

from .Vec import Vec


class Tensor:
    """Class that handles rank-2 tensor algebra."""

    def __init__(self, T):
        """
          Construct a rank-2 tensor
          Parameters
          ----------
            T : list
              list of four components of the tensor in the format Txx, Txy, Tyx, Tyy
        """
        if isinstance(T, list):
            self.T = np.array(T, dtype=np.float64)
        elif isinstance(T, np.ndarray):
            self.T = T 
        elif isinstance(T, Vec):
            TT = T.outer(T)
            self.T = np.array(TT.flatten(), dtype=np.float64)
        elif isinstance(T, Tensor):
            self.T = T.T 
        else:
            raise Exception('Unable to construct tensor.')

    def __add__(self, T):
        """
          Define Tensor addition.
          Parameters
          ----------
            T : Tensor
              Tensor to add to self
        """
        TT = self.T + T.T
        return Tensor(TT)

    def __sub__(self, T):
        """
          Define tensor subtraction.
          Parameters
          ----------
            T : Tensor
              Tensor to subtract from self
        """
        TT = self.T - T.T
        return Tensor(TT)

    def __mul__(self, scale):
        """
          Define Tensor scaling (from left).
          Parameters
          ----------
            scale : float
              Scaling factor
        """
        return Tensor(scale*self.T)

    def __rmul__(self, scale):
        """
          Define Tensor scaling (from right).
          Parameters
          ----------
            s : float
              Scaling factor
        """
        return Tensor(scale*self.T)

    def __iadd__(self, T):
        """
          Increment the Tensor by another Tensor.
          Parameters
          ----------
            T : Tensor
              Tensor to increment by
        """
        self.T += T.T
        return self

    def __isub__(self, T):
        """
          Decrement the Tensor by another vector.
          Parameters
          ----------
            T : Tensor
              Tensor to decrement by
        """
        self.T -= T.T
        return self

    def __repr__(self):
        """
          Return tensor components as a tuple
        """
        return "({:.6f}, {:.6f}, {:.6f}, {:.6f})".format(self.T[0], self.T[1], self.T[2], self.T[3])

    def __str__(self):
        """
          Return tensor components as a string for printing.
        """
        return "({:.6f}, {:.6f}, {:.6f}, {:.6f})".format(self.T[0], self.T[1], self.T[2], self.T[3])

    def dot(self, T):
        """
          Compute dot product between two Tensors (A : B, i.e. the double contraction).
          Parameters
          ----------
            T : Tensor
              Tensor to dot product with
        """
        A = self.T.reshape((2, 2))
        B = T.T.reshape((2, 2))
        return float(np.tensordot(A, B))

    def det(self):
        """
            Compute the determinant.
        """
        A = self.T.reshape((2, 2))
        return np.linalg.det(A)

    def trace(self):
        """
            Compute trace.
        """
        A = self.T.reshape((2, 2))
        return np.trace(A)

    def traceless(self):
        """
            Return traceless part of the tensor.
        """
        tr = self.trace()
        Txx, Txy, Tyx, Tyy = self.T
        return Tensor([Txx-0.5*tr, Txy, Tyx, Tyy - 0.5*tr])

    def eigenvals(self):
        """
            Return eigenvalues.
        """
        A = self.T.reshape((2, 2))
        return np.linalg.eig(A)[0]

    def eigenvectors(self):
        """
            Return eigenvectors.
        """
        A = self.T.reshape((2, 2))
        evecs = np.linalg.eig(A)[1]
        return (Vec(evecs[:, 0]), Vec(evecs[:, 1]))

    def principal_evec(self):
        """
            Return the eigenvector corresponding to the larger of the two eigenvalues. 
        """
        evs = self.eigenvals()
        evecs = self.eigenvectors()
        if evs[1] > evs[0]:
            return evecs[1]
        else:
            return evecs[0]
