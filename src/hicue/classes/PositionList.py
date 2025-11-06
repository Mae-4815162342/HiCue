from .imports import *

class PositionList():
    position_columns = ["Name", "Chromosome", "Start", "End", "Strand", "Tracks"]
    position_lock = threading.Lock()
    
    def __init__(self, nb_pos, high_low, binning = 1000):
        self._positions = {i:{} for i in range(nb_pos)}
        self._values = [np.inf if high_low == 'low' else -np.inf] * nb_pos
        self._high_low = high_low
        
    @staticmethod 
    def add_value_in_sorted_up(val, val_list, position, position_dict):
        """Adds a value in a ascending sorted list and returns the insertion index, or -1 if not inserted.
        Re-arranges the associated dictionnary."""
        if val <= val_list[0]:
            return -1
        val_list[0] = val
        position_dict[0] = position
        for i in range(len(val_list) - 1):
            if val_list[i] > val_list[i + 1]:
                val_list[i], val_list[i + 1] = val_list[i + 1], val_list[i]
                position_dict[i], position_dict[i + 1] = position_dict[i + 1], position_dict[i]
            else:
                break

    @staticmethod
    def add_value_in_sorted_down(val, val_list, position, position_dict):
        """Adds a value in a ascending sorted list and returns the insertion index, or -1 if not inserted.
        Re-arranges the associated dictionnary."""
        if val >= val_list[0]:
            return -1
        val_list[0] = val
        position_dict[0] = position
        for i in range(len(val_list) - 1):
            if val_list[i] <= val_list[i + 1]:
                val_list[i], val_list[i + 1] = val_list[i + 1], val_list[i]
                position_dict[i], position_dict[i + 1] = position_dict[i + 1], position_dict[i]
            else:
                break
        
    def append(self, position, value):
        """Appends a position to the position list at its value indexation in the values list."""
        with self.position_lock:
            match self._high_low:
                case "high":
                    value_tmp = value if value is not None else -np.inf
                    self.add_value_in_sorted_up(value_tmp, self._values, position, self._positions)
                case "low":
                    value_tmp = value if value is not None else np.inf
                    self.add_value_in_sorted_down(value_tmp, self._values, position, self._positions)
                
    def get_positions(self):
        """Returns the position list as a Panda dataframe."""
        df = pd.DataFrame.from_dict(self._positions, orient = "index")
        if "Index" in df.columns:
            df = df.set_index('Index')
        return df