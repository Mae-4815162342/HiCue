from hicue.classes.Reader import *
from hicue.classes.PairFormater import *
from hicue.classes.RandomSelector import *
from hicue.workers.MatrixExtractor import *

import time

test_files = [
    "SRP_dd_genes.bed", # bed file
    "TYs.gff", #gff file
    "plants_loops.bed2d" # bed2d file
]

expected_lengths = {
    "SRP_dd_genes.bed": 226, # track file
    "TYs.gff":53, #gff file
    "plants_loops.bed2d": 4 # both file
}

cool_files = {
    "SRP_dd_genes.bed": "Control.mcool", 
    "TYs.gff":"S_cerevisiae_TY.mcool",
    "plants_loops.bed2d": "Control.mcool"
}

format_parameters = {
    "SRP_dd_genes.bed": {
        "detrending": "ps",
        "circulars": ["NC_014500.1"],
        "diag_mask": 100,
    }, 
    "TYs.gff":{
        "detrending": "ps",
        "has_trans": True,
        "separate_by" : "chroms"
    },
    "plants_loops.bed2d": {
        "detrending": "patch",
        "circulars": ["NC_014500.1"],
        "diag_mask": 15000,
    }
}

test_parameters = {
    "SRP_dd_genes.bed": {
        "method":"mean", 
        "loop": False,
        "flip": True,
        "window":100000
    }, 
    "TYs.gff":{
        "method":"median", 
        "loop": True,
        "flip": False,
        "window":50000
    },
    "plants_loops.bed2d": {
        "method":"sum", 
        "loop": False,
        "flip": False,
        "window":150000
    }
}

nb_rand_per_pos = 2
window = 100000
binning = 5000
position_path = "test_data/positions/"
matrices_path = "test_data/matrices/"

### MatrixExtractor tests
for test_file in test_files:
    reader = Reader(position_path + test_file, test_file.split('.')[-1], loop = test_parameters[test_file]["loop"])
    positions, pair_queue = reader.read_file(threads=8)

    formater = PairFormater(positions, **format_parameters[test_file])
    formated_pairs = formater.format_pairs(pair_queue, threads = 8)

    cool_file = cooler.Cooler(f"{matrices_path}{cool_files[test_file]}::resolutions/{binning}")

    extraction_parameters = {
        "method":test_parameters[test_file]["method"],
        "nb_rand_per_pos":nb_rand_per_pos,
        "windows": [test_parameters[test_file]["window"]],
        "flip": test_parameters[test_file]["flip"]
    }

    start = time.time()
    pileup_queue = Queue()
    extractor = MatrixExtractor(formated_pairs, positions, **extraction_parameters)
    pileups = extractor.extract_from(cool_file, threads = 8)
    for pileup in pileups:
        pileup_queue.put(pileup)
    pileup_queue.put("DONE")
    stop = time.time()

    print(f"{test_file} extracted in {stop - start} seconds ({(stop - start) / 60} minutes).")

    random_queue = None
    if format_parameters[test_file]["detrending"] == "patch":
        extraction_parameters["randoms"] = True
        selector = RandomSelector(nb_rand_per_pos = nb_rand_per_pos)
        random_selection = selector.select_randoms(positions)

        random_queue = Queue()
        extractor = MatrixExtractor(random_queue, formated_pairs, random_selection, **extraction_parameters)
        pileups = extractor.extract_from(cool_file, randoms = False)
        for pileup in pileups:
            random_queue.put(pileup)
        random_queue.put("DONE")

    total_size = 0
    detrending_size = 0
    while True:
        try:
            pileup = pileup_queue.get()
            if random_queue != None:
                detrending = random_queue.get()
        except Empty:
            break
        if pileup == "DONE":
            break
        
        pileup_mat = pileup.get_matrix(test_parameters[test_file]["window"])

        if random_queue != None:
            detrending_mat = detrending.get_matrix(test_parameters[test_file]["window"])
            detrending_size += detrending.get_nb_matrices(test_parameters[test_file]["window"])

        total_size += pileup.get_nb_matrices(test_parameters[test_file]["window"])
    
    assert(total_size == len(formated_pairs))
    if random_queue != None:
        assert(detrending_size == len(formated_pairs) * nb_rand_per_pos)