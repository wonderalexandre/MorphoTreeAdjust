# MorphoTreeAdjust
MorphoTreeAdjust is a C++/Python implementation of a solution for the following problem:

Let $f$ be an image, and $\mathcal{T}^{\min}_{f}$ and $\mathcal{T}^{\max}_{f}$ be the min-tree and the max-tree of $f$ (built using the same connectivity).  Now, consider $L^{\max}\in\mathcal{T}^{\max}_{f}$ as a leaf in the max-tree $\mathcal{T}^{\max}_{f}$. Suppose $g=\Rec(\mathcal{T}_{f}^{\max} \setminus \{L^{\max}\})$ is the filtered image and consider $\mathcal{T}^{\min}_{g}$ and $\mathcal{T}^{\max}_{g}$ be the min-tree and the max-tree of $g$. 

MorphoTreeAdjust implements a solution for the following problem:
> What modifications do we need to make in the min-tree $\mathcal{T}^{\min}_{f}$ in order for the resulting min-tree to represent the same image as the max-tree $\mathcal{T}^{\max}_{g}$ ?