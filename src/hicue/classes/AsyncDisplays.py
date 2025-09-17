from .imports import *


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
        self._current_batch = []
        self._statics = static_params
        self._params_to_batch = params_to_batch
        self._nb_batch = 1
        self.start()
        
    def add_to_batch(self, **params):
        """Adds new object params to the current batch"""
        for param in self._params_to_batch:
            if param not in self._statics:
                self._statics[param] = params[param]
        self._current_batch.append(params)
        
    def empty_batch(self):
        self._current_batch = []
        
    def get_current_batch_size(self):
        return len(self._current_batch)
    
    def is_empty_batch(self):
        return self.get_current_batch_size() == 0
        
    def call_batch_display(self):
        asyncio.run(self._func(self._current_batch, batch_size = self._batch_size, title = f"batch#{self._nb_batch}", **self._statics))
        
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

            self.add_to_batch(**value)
            if self.get_current_batch_size() == self._batch_size:
                self.call_batch_display()
                self.empty_batch()
                self._nb_batch += 1
        if not self.is_empty_batch():
            self.call_batch_display()