import pandas as pd
import pickle
import numpy as np
import matplotlib.pyplot as plt
import gudhi


biodata=pd.read_csv('/Users/kave/Downloads/GSE200981_scRNAseq_processed.tsv', sep='\t', header=0)
with open('MDS.pickle', 'rb') as f:
    bioCoord = pickle.load(f)
print(np.mean(bioCoord,axis=0))


# Sample data (replace with your actual data)
# Create the scatter plot
plt.scatter(x=bioCoord[:,0],  y=bioCoord[:,1])  # 's' controls marker size
plt.grid(True)

# Show the plot
plt.show()

rips = gudhi.RipsComplex(points=bioCoord, max_edge_length=42)
simplex_tree = rips.create_simplex_tree(max_dimension=1)


diag = simplex_tree.persistence(homology_coeff_field=2, min_persistence=0)
print("diag=", diag)

pplot = gudhi.plot_persistence_diagram(diag)
pplot.show()

