from hicue.classes.Reader import *
import time

# reader test
performance_tests = [
    "test_data/positions/some_human_genes.bed", 
    #"test_data/gff/human_hg38.gff"
] # test for the maximum input size (to be tested with the tracks selection later

for file in performance_tests:
    start = time.time()
    reader = Reader(file, file.split('.')[-1], record_type = "gene")
    pos, _ = reader.read_file(threads=8)
    end = time.time()
    print(file, end - start)
    print(pos.head())