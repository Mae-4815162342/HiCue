import cooler
import matplotlib.pyplot as plt
import numpy as np

matrix = "/data/Maelys/datas/cool_files/Control/Control.mcool::resolutions/5000"
cool = cooler.Cooler(matrix)

plt.imshow(np.log10(cool.matrix(balance=True)[:]), extent=[0, cool.chromsizes[cool.chromnames[0]], cool.chromsizes[cool.chromnames[0]], 0], cmap = "afmhot_r")
plt.show()