# MorphoTreeAdjust
MorphoTreeAdjust is a C++/Python implementation of a solution for the following problem:

Let $f$ be an image, and $T^\min_f$ and $T^\max_f$ be the min-tree and the max-tree of $f$.
We set:
- $L$, a leaf node in the min-tree $T^\min_f$;
- $T^\min_g$, the min-tree obtained after pruning the leaf node $L$;
- $g$, the image reconstructed from the pruned min-tree $T^\min_g$;
- $T^\max_g$ the max-tree of the image $g$.  

Problem:
> How to modify $T^\max_f$ to obtain $T^\max_g$?

Demo:
> [Click here and see a demo on jupyter notebook](./notebooks/morphoTreeAdjust_example_leaf.ipynb)

Paper:
Alves, W. A. L., Passat, N., Silva, D. J., Morimitsu, A., & Hashimoto, R. F. (2025). *Efficient connected alternating sequential filters based on component trees*. Paper submitted to the 4th International Conference on Discrete Geometry and Mathematical Morphology (DGMM 2025).
