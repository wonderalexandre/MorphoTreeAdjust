import numpy as np

import morphoTreeAdjust as mta


image = np.array(
    [
        [4, 4, 2, 1],
        [4, 3, 2, 1],
        [5, 5, 2, 0],
        [5, 6, 6, 0],
    ],
    dtype=np.uint8,
)

adj = mta.AdjacencyRelation(image.shape[0], image.shape[1], 1.5)
maxtree = mta.DynamicComponentTree(image, True, adj)
mintree = mta.DynamicComponentTree(image, False, adj)

adjust = mta.DynamicComponentTreeAdjustment(mintree, maxtree)
adjust.refreshAttributes()

casf = mta.ComponentTreeCasf(image, "area", adj)
filtered = casf.filter([1, 2])

print("num nodes (max-tree):", maxtree.numNodes)
print("filtered image:")
print(filtered)
