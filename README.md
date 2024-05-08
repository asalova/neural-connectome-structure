# neural-connectome-structure
Repository for the manuscript "Combined topological and spatial constraints are required to capture the structure of neural connectomes"

The folders titled *_data (where * is fly, mouse, and human) contain the processed individual neuron information and connectome/contactome structure.
Below is a brief summary of the columns of the .csv files contained in each of these folders.

# Summary of processed data

## Neuron information

### *_basic_neuron_info.csv
**id** : ids of the neurons in the publicly released data sets <br />
**x_cm (y_cm, z_cm)** : x (y, z) position of the center of mesh, nm

### *_extended_neuron_info.csv
**id** : ids of the neurons in the publicly released data sets <br />
**x_soma (y_soma, z_soma)** : x (y, z) position of the soma, nm <br />
**x_cm (y_cm, z_cm)** : x (y, z) position of the center of mesh, nm <br />
**degree** : undirected degree of each neuron <br />
**presynaptic_degree** : pre-synaptic degree (# of post-synaptic neighbors) of each neuron <br />
**postsynaptic_degree** : post-synaptic degree (# of post-synaptic neighbors) of each neuron <br />
**weighted_presynaptic_degree** : # of relevant (established with other neurons we consider) pre-synapses of each neuron <br />
**weighted_postsynaptic_degree** : # of relevant (established with other neurons we consider) post-synapses of each neuron <br />
**pca_1_x (pca_1_y, pca_1_z)** : dominant principle component x (y, z) direction based on the mesh vertices <br />
**evr_1** : explained variance ratio of the dominant principle component direction <br />
**linear_span_pca_1** : linear span of the neuron along the dominant principle component direction  <br />
**n_mesh_vertices** : number of mesh vertices

Note that we use the index of the row corresponding to each neuron (e.g., integers from 0 to 15731 for the human data) as connectome and contactome node labels in the .csv files below.

## Connectome

### *_directed_connectome.csv 
Here, each row corresponds to a pair of neurons **(i,j)**, s.t. an edge from i to j exists <br />
**i** : pre-synaptic neuron index <br />
**j** : post-synaptic neuron index <br />
**weight** : number of synapses from neuron i to neuron j

### *_weighted_connectome.csv
Here, each row corresponds to a pair of neurons **(i,j)**, i<j, s.t. an edge between i and j exists <br />
**i** : neuron index <br />
**j** : neuron index <br />
**reciprocated** : indicates whether the directed edge between i and j is reciprocated (exists in both directions)

### *_connectome_edge_weights.csv
**weight** : number of synapses between two neurons <br />
**count** : number of pairs of neurons with an edge of a given weight

## Contactome

### *_connectome_edge_weights.csv
Here, each row corresponds to a pair of neurons in physical contact <br />
**i** : neuron index <br />
**j** : neuron index <br />

# Summary of models
We include the edge probabilities obtained for each model discussed in the manuscript (ER, model_c, model_d, model_d_c, model_k, model_k_c, model k_L) as .npy files in models/*.zip, where * is fly, mouse, and human.  <br />
As an example of getting edge probabilities in Python: <br />
```python
import numpy as np
from scipy.spatial.distance import squareform

p = squareform(np.load('models/mouse/p_model_k.npy'))
```
There, 'p[i,j]' correpo
