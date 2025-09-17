from hicue.classes.Pileup import Pileup
import numpy as np
import time
import matplotlib.pyplot as plt

# testing pileup capacity

for k in range(1000, 46000, 1000):
    pileup = Pileup(nb_matrices=k)

    start = time.time()
    for i in range(k):
        if i % 1000 == 0:
            print(i)
        matrix = np.random.uniform(low=0, high=10, size=(201, 201))
        np.fill_diagonal(matrix, i)
        pileup.add_submatrix(100000, matrix)
    print("It took ", time.time() - start, " seconds append ", k, " matrices to the pileup.")

    start = time.time()
    pileup_matrix = pileup.get_matrix(100000)
    print(np.nanmedian(range(k)), np.diag(pileup_matrix)[0])
    print("It took ", time.time() - start, " seconds to retrieve the pileup of ", k, " matrices.")
