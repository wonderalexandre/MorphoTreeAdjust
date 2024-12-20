# MorphoTreeAdjust
MorphoTreeAdjust is a C++/Python implementation of a solution for the following problem:

Let $f$ be an image, and $T^\min_f$ and $T^\max_f$ be the min-tree and the max-tree of $f$ (built using the same connectivity).  Now, consider $L^\max \in T^\max_f$ as a leaf in the max-tree $T^\max_f$. Suppose $g=\text{Rec}(T^\max_f \setminus \{L^\max\})$ is the filtered image and consider $T^\min_g$ and $T^\max_g$ be the min-tree and the max-tree of $g$. 

MorphoTreeAdjust implements a solution for the following problem:
> What modifications do we need to make in the min-tree $T^\min_f$ in order for the resulting min-tree to represent the same image as the max-tree $T^\max_g$ ?
