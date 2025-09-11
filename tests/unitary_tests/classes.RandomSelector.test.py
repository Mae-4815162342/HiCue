from hicue.classes.Reader import *
from hicue.classes.RandomSelector import *

# TODO: re-write execution later, for now test with :
# ../.venv/bin/python3.10 /data/Maelys/visualisation_tools/hicue/tests/unitary_tests/classes.RandomSelector.test.py

## Random positions selection test.

test_files = [
    "SRP_dd_genes.bed", # bed file
    "TYs.gff", #gff file
    "plants_loops.bed2d" # bed2d file
]
position_path = "test_data/positions/"
nb_rand_per_pos = 2

for test_file in test_files:
    reader = Reader(position_path + test_file, test_file.split('.')[-1])
    positions, pair_queue = reader.read_file(threads=8)

    selector = RandomSelector(nb_rand_per_pos = nb_rand_per_pos)
    random_selection = selector.select_randoms(positions)

    random_selection.to_csv(f"test_out/random_pos_{test_file.split('.')[0]}")

    assert(len(random_selection) == len(positions) * nb_rand_per_pos)
    assert(len(random_selection) == len(np.unique(random_selection.index)))
