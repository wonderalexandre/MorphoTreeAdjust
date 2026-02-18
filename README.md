# MorphoTreeAdjust
MorphoTreeAdjust is a C++/Python implementation of a solution for the following problem:

Let $f$ be an image, and $T^\min_f$ and $T^\max_f$ be the min-tree and the max-tree of $f$.
We set:
- $S$, a subtree in the min-tree $T^\min_f$;
- $T^\min_g$, the min-tree obtained after pruning the subtree $S$ from $T^\min_f$;
- $g$, the image reconstructed from the pruned min-tree $T^\min_g$;
- $T^\max_g$ the max-tree of the image $g$.  

Problem:
> How to modify $T^\max_f$ to obtain $T^\max_g$?

<br>
### Demo:
> [Click here and see a demo on jupyter notebook - leaf](./notebooks/morphoTreeAdjust_example_leaf.ipynb) <br>
> [Click here and see a demo on jupyter notebook - subtree](./notebooks/morphoTreeAdjust_example_subtree.ipynb)

### Paper:
> Wonder Alves, Nicolas Passat, Dênnis José da Silva, Alexandre Morimitsu, Ronaldo F. Hashimoto. *Efficient connected alternating sequential filters based on component trees*. International Conference on Discrete Geometry and Mathematical Morphology (DGMM), Nov 2025, Groningen, Netherlands. ⟨[hal-05163556](https://hal.science/hal-05163556/)⟩

> [Under review] Wonder Alves, Nicolas Passat, Dênnis José da Silva, Alexandre Morimitsu, Ronaldo F. Hashimoto. *Component tree: Update rather than rebuild*. Journal of Mathematical Imaging and Vision, 2026.


### Benchmark
Build and run the C++ benchmarks:

```sh
cd tests
cmake -S . -B build
cmake --build build
./build/jmiv2026_benchmark --repeat 5 --warmup 1 --json --no-validate image1.png image2.png
```
