from .imports import *

def schedule_workers(worker_class, worker_location, threads, **wargs):
    """Creates a worker of class worker_class per thread with the arguments wargs. Returns the instances created in a list."""
    worker_module = importlib.import_module(worker_location)
    worker_module.initialize_globals()
    WorkerClass = getattr(importlib.import_module(worker_location), worker_class)
    workers = []
    for _ in range(threads):
        workers.append(WorkerClass(**wargs))
    return workers

def join_workers(workers):
    """Joins workers in the workers list."""
    for worker in workers:
        worker.join()
        
def join_queues(queues, threads = 1):
    """Signal DONE in each queue of queues"""
    for queue in queues:
        for _ in range(threads):
            queue.put("DONE")