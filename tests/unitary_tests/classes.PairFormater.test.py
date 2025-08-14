from hicue.classes.Reader import *
from hicue.classes.PairFormater import *


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

position_path = "test_data/positions/"

annotations = {
    "SRP_dd_genes.bed": {'tracks':"test_data/tracks/PGA_stat.bw"}, # track file
    "TYs.gff":{'gff':"test_data/gff/W303.gff3"}, #gff file
    "plants_loops.bed2d": {'tracks':"test_data/tracks/PGA_stat.bw",'gff':"test_data/gff/D_dadantii.gff"} # both file
}

### PairFormater tests

test_file = "TYs.gff" # using multichromosomal positions in loop mode for testing trans/cis contacts

reader = Reader(position_path + test_file, test_file.split('.')[-1], annotation_files=annotations[test_file], loop = True)
positions, pair_queue = reader.read_file(threads=8)

# building one pair queue per test
nb_test = 7
queues = [Queue() for _ in range(nb_test)]
nb_pairs = 0
nb_cis = 0
while True:
    try: 
        val = pair_queue.get(timeout=2)
        if val != "DONE":
            nb_pairs += 1
            if positions.loc[val[0]]["Chromosome"] == positions.loc[val[1]]["Chromosome"]:
                nb_cis += 1
        for queue in queues:
            queue.put(val)
    except Empty:
        break

# no parameters
formater = PairFormater(positions)
formated_pairs = formater.format_pairs(queues[0], threads = 8)
assert(len(formated_pairs) == nb_cis)

# static parameters
test_params = {
    "center" : "start",
    "detrending": "ps",
    "diag_mask": 10000,
    "has_trans": True,
    "circulars": ["chrI"]
}
formater = PairFormater(positions, **test_params)
formated_pairs = formater.format_pairs(queues[1], threads = 8)
assert(len(formated_pairs) == nb_pairs)
for _, pair in formated_pairs.head(10).iterrows(): # checking circularity
    chrom1 = positions.loc[pair["Locus1"]]["Chromosome"]
    chrom2 = positions.loc[pair["Locus2"]]["Chromosome"]
    assert((chrom1 == "chrI") == pair["Chrom1_circular"])
    assert((chrom2 == "chrI") == pair["Chrom2_circular"])

# separate by chromosomes
test_params = {
    "separate_by" : "chroms",
    "center" : "start",
    "detrending": "ps",
    "diag_mask": 10000,
    "has_trans": True
}
formater = PairFormater(positions, **test_params)
formated_pairs = formater.format_pairs(queues[2], threads = 8)
assert(len(formated_pairs) == nb_pairs)

# separate by direction
test_params = {
    "separate_by" : "chroms,direct",
    "center" : "start",
    "detrending": "ps",
    "diag_mask": 10000,
    "has_trans": True
}
formater = PairFormater(positions, **test_params)
formated_pairs = formater.format_pairs(queues[3], threads = 8)
assert(len(formated_pairs) == nb_cis) # only cis contact can be separated by direction

# separate by regions
test_params = {
    "separate_by" : "regions",
    "center" : "start",
    "separate_regions" : "test_data/others/regions_test.csv",
    "detrending": "ps",
    "diag_mask": 10000,
    "has_trans": True
}
formater = PairFormater(positions, **test_params)
formated_pairs = formater.format_pairs(queues[4], threads = 8)
assert(len(formated_pairs) == 3) # only three pairs in the dataset are contained in the delimited regions as pairs

# separate by contact_range
test_params = {
    "center" : "start",
    "contact_range" : [200000, 10000000, 100000],
    "detrending": "ps",
    "diag_mask": 10000,
    "has_trans": True
}
formater = PairFormater(positions, **test_params)
formated_pairs = formater.format_pairs(queues[5], threads = 8)
for _, pair in formated_pairs.head(10).iterrows():
    interval = pair["Sep_id"].split("_")[0]
    start, stop = int(interval.split("-")[0]), int(interval.split("-")[1])
    assert(pair["Distance"] >= start)
    assert(pair["Distance"] < stop)

# separate by all + min_dist
test_params = {
    "separate_by" : "chroms,direct,regions",
    "center" : "start",
    "contact_range" : [200000, 10000000, 100000],
    "separate_regions" : "test_data/others/regions_test.csv",
    "min_dist" : 15000,
    "detrending": "ps",
    "diag_mask": 10000,
    "has_trans": True
}
formater = PairFormater(positions, **test_params)
formated_pairs = formater.format_pairs(queues[6], threads = 8)
assert(len(formated_pairs) == 1) # only one pair fills all conditions