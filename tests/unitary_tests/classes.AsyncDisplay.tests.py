from hicue.classes.AsyncDisplays import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if not os.path.exists("test_out/displays"):
    os.mkdir("test_out/displays")

## Display class test
# parameters
matrix = np.random.random((100,100))
display_queue = Queue()
nb_samples = 4

# test method
async def display_matrix(matrix, matrix_name, cmap = "bwr", outpath = ""):
    plt.figure()
    plt.imshow(matrix, cmap = cmap)
    plt.title(matrix_name)
    plt.savefig(f"{outpath}/{matrix_name}.png")
    
# static parameters
displayer = Display(
    input_queue = display_queue,
    output_queues=[],
    function = display_matrix,
    cmap = "afmhot_r",
    outpath = "test_out/displays"
)

for i in range(nb_samples):
    display_queue.put({
        "matrix": matrix, 
        "matrix_name": f"test_{ i + 1}"
    })
display_queue.put("DONE")

# assertions
displayer.join()
for i in range(nb_samples):
    assert(os.path.exists(f"test_out/displays/test_{i + 1}.png"))

## DislayBatch class test
# parameters
matrix = np.random.random((100,100))
display_queue = Queue()
nb_samples = 15
batch_size = 8

# test method
async def display_matrices(matrices, batch_size = 8, title = "", cmap = "bwr", outpath = ""):
    plt.figure()    
    nb_matrices = len(matrices)
    
    cols = math.ceil(math.sqrt(nb_matrices))
    rows = math.ceil(nb_matrices / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 3*rows))
    axes = np.array(axes).reshape(-1)
    for i in range(nb_matrices):
        axes[i].imshow(matrices[i]["matrix"], cmap = cmap)
        axes[i].set_title(f'{matrices[i]["matrix_name"]}')

    # Supprimer les subplots vides (si nb_matrices < rows * cols)
    for j in range(nb_matrices, rows * cols):
        fig.delaxes(axes[j])

    fig.tight_layout()
    plt.savefig(f"{outpath}/{title}.png")
    
# static parameters
displayer = DisplayBatch(
    input_queue = display_queue,
    output_queues=[],
    function = display_matrices,
    batch_size = batch_size,
    cmap = "afmhot_r",
    outpath = "test_out/displays"
)

if not os.path.exists("test_out/displays"):
    os.mkdir("test_out/displays")

for i in range(nb_samples):
    display_queue.put({
        "matrix": matrix, 
        "matrix_name": f"test_{ i + 1}"
    })
display_queue.put("DONE")

displayer.join()
for i in range(1, displayer._nb_batch):
    assert(os.path.exists(f"test_out/displays/batch#{i}.png"))