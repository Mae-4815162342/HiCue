from hicue.classes.Pileup import *

# initialisation test
pileup = Pileup()
assert(len(np.unique(pileup.get_matrix(11))) == 1)

# update test
pileup.add_value_to_pixel(8, 9, 10.5)
pileup.add_value_to_pixel(4, 6, 5.5)
pileup.add_value_to_pixel(0, 8, 3.5)
submatrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
pileup.add_submatrix(submatrix)

assert(np.nansum(pileup.get_matrix(11)) == 10.5  + 5.5 + 3.5 + np.sum(submatrix))

# multithreading test
pileup = Pileup(mode = "median")
rand_list = [random.random() for i in range(1000)]
workers = []
for value in rand_list:
    t = threading.Thread(target = pileup.add_submatrix, args = (np.array([[value]]),))
    t.start()
    workers.append(t)
for w in workers:
    w.join()

assert(np.nanmedian(rand_list) == pileup.get_matrix(1)[0, 0])

pileup = Pileup(mode = "mean")
rand_list = [random.random() for i in range(1000)]
workers = []
for value in rand_list:
    t = threading.Thread(target = pileup.add_submatrix, args = (np.array([[value]]),))
    t.start()
    workers.append(t)
for w in workers:
    w.join()

assert(np.nanmean(rand_list) == pileup.get_matrix(1)[0, 0])

pileup = Pileup(mode = "sum")
rand_list = [random.random() for i in range(1000)]
workers = []
for value in rand_list:
    t = threading.Thread(target = pileup.add_submatrix, args = (np.array([[value]]),))
    t.start()
    workers.append(t)
for w in workers:
    w.join()

assert(np.nansum(rand_list) == pileup.get_matrix(1)[0, 0])