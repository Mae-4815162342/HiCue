from .imports import *

class Pileup():    
    pileup_lock = threading.Lock()
    
    def __init__(self, nb_matrices = 0, sep_id = "", mode = "median", binning = 1000, cool_name = ""):
        self._max_matrices = nb_matrices
        self._sep_id = sep_id
        self._cool_name = cool_name
        self._mode = mode
        self._binning = binning
        self._pileup_matrices = {}
        self._nb_matrices = {}
        self._size = {}

    def add_submatrix(self, window, submatrix):
        """Adds a submatrix to the pileup."""
        sub_size = len(submatrix)
        with self.pileup_lock:
            if window not in self._pileup_matrices:
                self._pileup_matrices[window] = np.empty((self._max_matrices, sub_size**2), dtype=np.float32)
                self._size[window] = sub_size
                self._nb_matrices[window] = 0
            flattened = submatrix.flatten()
            self._pileup_matrices[window][self._nb_matrices[window]] = flattened
            self._nb_matrices[window] += 1

    def get_matrix(self, window):
        """Returns the numpy array of the pileup of matrix_size."""
        if self._size[window] > 0:
            matrix = None
            match self._mode:
                case "median":
                    matrix = np.nanmedian(self._pileup_matrices[window], axis = 0)
                case "mean":
                    matrix = np.nanmean(self._pileup_matrices[window], axis = 0)
                case "sum":
                    matrix = np.nansum(self._pileup_matrices[window], axis = 0)
            return matrix.reshape((self._size[window], self._size[window]))
        return None
    
    def get_sep_id(self):
        """Getter for sep_id"""
        return self._sep_id
    
    def get_cool_name(self):
        """Getter for cool_name"""
        return self._cool_name
    
    def get_binning(self):
        """Getter for binning"""
        return self._binning
    
    def get_size(self, window):
        """Getter for size"""
        return self._size[window]
    
    def get_nb_matrices(self, window):
        """Getter for nb_matrices"""
        nb_matrices = 0
        with self.pileup_lock:
            nb_matrices = self._nb_matrices[window]
        return nb_matrices