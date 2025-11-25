from hicue.utils import *

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
        self._tracks = {}
        self._directions = {}

    def add_submatrix(self, window, submatrix, directions = [1, 1]):
        """Adds a submatrix to the pileup."""
        sub_size = len(submatrix[0])
        tracks = len(submatrix) - sub_size
        with self.pileup_lock:
            if window not in self._pileup_matrices:
                self._pileup_matrices[window] = np.empty((self._max_matrices, sub_size**2 + sub_size * tracks), dtype=np.float32)
                self._size[window] = sub_size
                self._tracks[window] = tracks
                self._nb_matrices[window] = 0
                self._directions[window] = directions
            self._pileup_matrices[window][self._nb_matrices[window]] = submatrix.flatten()
            self._nb_matrices[window] += 1

    def get_matrix(self, window):
        """Returns the numpy array of the pileup of matrix_size."""
        with self.pileup_lock:
            if window not in self._size:
                return None
            if self._size[window] > 0:
                matrix = None
                match self._mode:
                    case "median":
                        matrix = np.nanmedian(self._pileup_matrices[window][:self._nb_matrices[window]], axis = 0)
                    case "mean":
                        matrix = np.nanmean(self._pileup_matrices[window][:self._nb_matrices[window]], axis = 0)
                    case "sum":
                        matrix = np.nansum(self._pileup_matrices[window][:self._nb_matrices[window]], axis = 0)
                return matrix.reshape((self._size[window] + self._tracks[window], self._size[window]))
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

    def has_tracks(self, window):
        """Getter for size"""
        return self._tracks[window] > 0
    
    def is_directed(axis=None):
        """Returns the direction existence depending of the axis if provided."""
    
    def get_nb_matrices(self, window):
        """Getter for nb_matrices"""
        nb_matrices = 0
        with self.pileup_lock:
            nb_matrices = self._nb_matrices[window]
        return nb_matrices