# MorphoTreeAdjust
MorphoTreeAdjust is a C++/Python implementation of a solution for the following problem:

Let $f$ be an image, and $T^{\min}_f$ and $T^{\max}_f$ be the min-tree and the max-tree of $f$ (built using the same connectivity).  Now, consider $L^{\max} \in T^{\max}_{f}$ as a leaf in the max-tree $T^{\max}_{f}$. Suppose $g=\text{Rec}(T^{\max}_{f} \setminus \{L^{\max}\})$ is the filtered image and consider $T^{\min}_{g}$ and $T^{\max}_{g}$ be the min-tree and the max-tree of $g$. 

MorphoTreeAdjust implements a solution for the following problem:
> What modifications do we need to make in the min-tree $T^{\min}_{f}$ in order for the resulting min-tree to represent the same image as the max-tree $T^{\max}_{g}$ ?
