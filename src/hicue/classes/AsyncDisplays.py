from hicue.imports import *


def initialize_globals():
    pass

# the asynchronous displayer require matplotlib.use('Agg') (before any other import)

class Display(threading.Thread):
    def __init__(self, input_queue, output_queues, function, **static_params):
        super(Display, self).__init__()
        self._input_queue = input_queue 
        self._output_queues = output_queues
        self._func = function
        self._statics = static_params
        self.start()
        
    def call_display(self, **params):
        asyncio.run(self._func(**params, **self._statics))
        
    def run(self):
        while True:
            try: 
                value = self._input_queue.get()
                for queue in self._output_queues:
                    queue.put(value)
            except Empty:
                break
            if value == "DONE":
                break
            params = value
            self.call_display(**params)

class DisplayBatch(threading.Thread):
    """Same display calls but will batch the inputs with batch size."""
    def __init__(self, input_queue, output_queues, function, batch_size, params_to_batch = [], **static_params):
        super(DisplayBatch, self).__init__()
        self._input_queue = input_queue
        self._output_queues = output_queues
        self._func = function
        self._batch_size = batch_size
        self._current_batch = {}
        self._statics = static_params
        self._params_to_batch = params_to_batch
        self._outfolders = {}
        self._keys = []
        self._nb_batch = {}
        self.start()
        
    def add_to_batch(self, window, sep_id, **params):
        """Adds new object params to the current batch"""
        for param in self._params_to_batch:
            if param not in self._statics:
                self._statics[param] = params[param]
        if f"{window}.{sep_id}" not in self._current_batch:
            self._current_batch[f"{window}.{sep_id}"] = []
            self._nb_batch[f"{window}.{sep_id}"] = 1
            self._outfolders[f"{window}.{sep_id}"] = params['outfolder']
            self._keys.append((window,sep_id))
        self._current_batch[f"{window}.{sep_id}"].append(params)
        
    def empty_batch(self, window, sep_id):
        self._current_batch[f"{window}.{sep_id}"] = []
        
    def get_current_batch_size(self, window, sep_id):
        return len(self._current_batch[f"{window}.{sep_id}"])
    
    def is_empty_batch(self, window, sep_id):
        return self.get_current_batch_size(window, sep_id) == 0
        
    def call_batch_display(self, window, sep_id):
        asyncio.run(self._func(self._current_batch[f"{window}.{sep_id}"], batch_size = self._batch_size, outfolder = self._outfolders[f"{window}.{sep_id}"], title = f"batch#{self._nb_batch[f'{window}.{sep_id}']}", window = window, **self._statics))
        
    def run(self):
        while True:
            try: 
                value = self._input_queue.get()
                for queue in self._output_queues:
                    queue.put(value)
            except Empty:
                break
            if value == "DONE":
                break

            window = value["window"]
            sep_id = value['pair']["Sep_id"]
            self.add_to_batch(**value, sep_id=sep_id)
            if self.get_current_batch_size(window, sep_id) == self._batch_size:
                self.call_batch_display(window, sep_id)
                self.empty_batch(window, sep_id)
                self._nb_batch[f"{window}.{sep_id}"] += 1
        for window, sep_id in self._keys:
            if not self.is_empty_batch(window, sep_id):
                self.call_batch_display(window, sep_id)