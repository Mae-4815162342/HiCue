from hicue.workers.utils import *

# testing extract_matrix
matrix_path = "test_data/matrices/Control.mcool"
chromosome = "NC_014500.1"
windows = [30000, 50000, 100000]
resolutions = [1000, 5000, 10000]

# 3 positions
# inferior overflow
inf = {
    "Chromosome":chromosome,
    "Start":4430,
    "End":6299,
    "Strand":"+"
}
# superior overflow
sup = {
    "Chromosome":chromosome,
    "Start":491000,
    "End":4918841,
    "Strand":"-"
}
# no overflow
center = {
    "Chromosome":chromosome,
    "Start":3243473,
    "End":3244109,
    "Strand":"-"
}

test_positions = [
    (center, center), # no overflow
    (inf, center), # dim 1 inf
    (sup, center), # dim 1 sup
    (center, inf), # dim 2 inf
    (center, sup), # dim 2 sup
    (inf, inf), # dim 1 inf/dim 2 inf
    (inf, sup), # dim 1 inf/dim 2 sup
    (sup, inf), # dim 1 sup/dim 2 inf
    (sup, sup)  # dim 1 sup/dim 2 sup
]

for res in resolutions:
    cool_file = cooler.Cooler(f"{matrix_path}::resolutions/{res}")
    for window in windows:
        
        expected_size = (window//res) * 2 + 1
        
        for pos1, pos2 in test_positions:
            sub = extract_window(cool_file, pos1, pos2, res, window)
            sub_center = extract_window(cool_file, pos1, pos2, res, window, center="center")
            sub_end = extract_window(cool_file, pos1, pos2, res, window, center="end")
            sub_raw = extract_window(cool_file, pos1, pos2, res, window, raw=True)
            sub_circ = extract_window(cool_file, pos1, pos2, res, window, is_loc1_circ= True, is_loc2_circ= True)
            
            assert(sub.shape[0] == expected_size)            
            assert(sub.shape[1] == expected_size)
            assert(sub_center.shape[0] == expected_size)            
            assert(sub_center.shape[1] == expected_size)
            assert(sub_end.shape[0] == expected_size)            
            assert(sub_end.shape[1] == expected_size)
            assert(sub_raw.shape[0] == expected_size)            
            assert(sub_raw.shape[1] == expected_size)
            assert(sub_circ.shape[0] == expected_size)            
            assert(sub_circ.shape[1] == expected_size)

            #fig, (ax1, ax2) = plt.subplots(1, 2)
            #fig.suptitle(f"{pos1}\n{pos2}")
            #ax1.matshow(sub)
            #ax1.set_title("Filled with NaN")
            #ax2.matshow(sub_circ)
            #ax2.set_title("Filled with circular values")
            #plt.show()

# testing compute_distance
new_position = {
    "Chromosome":chromosome,
    "Start":center["Start"] + 50000,
    "End":center["End"] + 50000,
    "Strand":"-"
}
new_position2 = {
    "Chromosome":chromosome,
    "Start":center["Start"] - 100000,
    "End":center["End"] - 100000,
    "Strand":"-"
}

dist1 = compute_distance(center, center)
dist2 = compute_distance(center, new_position)
dist3 = compute_distance(center, new_position2)

assert(dist1 == 0)
assert(dist2 == 50000)
assert(dist3 == - 100000)

# testing mask_diagonal
diag_sizes = [1000, 10000, 50000]
for res in resolutions:
    cool_file = cooler.Cooler(f"{matrix_path}::resolutions/{res}")
    for window in windows:
        for diagonal_mask in diag_sizes:
            sub1 = extract_window(cool_file, center, center, res, window)
            mask_diagonal(sub1, center, center, res, diagonal_mask)
            sum1 = 0
            for i in range(diagonal_mask//res):
                val1 = np.nansum(np.diag(sub1[i:]))
                val2 = np.nansum(np.diag(sub1[:, i:]))
                sum1 = np.nansum([val1, val2, sum1])
            assert(sum1 == 0)   
            

            sub2 = extract_window(cool_file, center, new_position, res, window)
            mask_diagonal(sub2, center, new_position, res, diagonal_mask)
            sum2 = 0
            for i in range(diagonal_mask//res):
                val1 = np.nansum(np.diag(sub2[abs(dist2) + i:, :-(abs(dist2) + i)]))
                val2 = np.nansum(np.diag(sub2[abs(dist2) - i:, :- (abs(dist2) - i)]))
                sum2 = np.nansum([val1, val2, sum2])
            assert(sum2 == 0)
            
            sub3 = extract_window(cool_file, center, new_position2, res, window)
            mask_diagonal(sub3, center, new_position2, res, diagonal_mask)
            sum3 = 0
            for i in range(diagonal_mask//res):
                val1 = np.nansum(np.diag(sub3[:- (abs(dist3) + i),  abs(dist3) + i:]))
                val2 = np.nansum(np.diag(sub3[:-  (abs(dist3) - i),  abs(dist3) - i:]))
                sum3 = np.nansum([val1, val2, sum3])
            assert(sum3 == 0)

# testing detrending_submatrix
expected_sums = [
    2917.3826767482274,
    3366.8404795326105,
    4181.701936240688,
    9486.869747035627,
    10041.12354776881,
    11183.366924694676,
    39590.888631871145,
    40659.798296686175,
    40751.25743580216,
    138.23066609539777,
    156.5361022066772,
    193.89610490103678,
    421.85490210772616,
    446.6459250206324,
    494.38218885796874,
    1667.1769196884952,
    1733.2602408281557,
    1725.2775044009895,
    40.03409779614635,
    46.93417935463238,
    54.348473575902744,
    116.41741207575299,
    123.14832387707187,
    135.09028693594288,
    436.7332914270096,
    453.52285277210757,
    451.4906370397247
]

k = 0
for res in resolutions:
    cool_file = cooler.Cooler(f"{matrix_path}::resolutions/{res}")
    ps = distance_law(cool_file.matrix().fetch(chromosome))
    ps[np.isnan(ps)] = 0.0
    
    for window in windows:
        sub1 = extract_window(cool_file, center, center, res, window)
        sub1_det = detrend_submatrix(sub1, center, center, res, ps)
        assert(abs(np.nansum(sub1_det)- expected_sums[k]) < 1e-9)
        k += 1
        
        #fig, (ax1, ax2) = plt.subplots(1, 2)
        #fig.suptitle(f"Detrending center-center")
        #ax1.matshow(sub1)
        #ax1.set_title("Normalized")
        #ax2.matshow(sub1_det)
        #ax2.set_title("Detrended")
        #plt.show()
        
        sub2 = extract_window(cool_file, center, new_position, res, window)
        sub2_det = detrend_submatrix(sub2, center, new_position, res, ps)
        assert(abs(np.nansum(sub2_det)- expected_sums[k]) < 1e-9)
        k += 1
        
        #fig, (ax1, ax2) = plt.subplots(1, 2)
        #fig.suptitle(f"Detrending center-new_position")
        #ax1.matshow(sub2)
        #ax1.set_title("Normalized")
        #ax2.matshow(sub2_det)
        #ax2.set_title("Detrended")
        #plt.show()
        
        sub3 = extract_window(cool_file, center, new_position2, res, window)
        sub3_det = detrend_submatrix(sub3, center, new_position2, res, ps)
        assert(abs(np.nansum(sub3_det)- expected_sums[k]) < 1e-9)
        k += 1
        
        #fig, (ax1, ax2) = plt.subplots(1, 2)
        #fig.suptitle(f"Detrending center-new_position2")
        #ax1.matshow(sub3)
        #ax1.set_title("Normalized")
        #ax2.matshow(sub3_det)
        #ax2.set_title("Detrended")
        #plt.show()