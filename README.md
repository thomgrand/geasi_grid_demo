# GEASI Demo
This repository contains a minimal working example of our algorithm presented in GEASI: Geodesics-based Earliest Activation Sequence
Identification in cardiac eikonal models *TODO: Link*. The main steps of the algorithm are similar, but to offer a concise and short code some simplifying assumptions were made (see Limitations).

The method is demonstrated in [GEASI_Grid.ipynb](GEASI_Grid.ipynb), whereas [utils.py](utils.py) contains some utility functions to keep the code short. The presented problem is similar to the square-domain example of the paper with slightly varying results due to the deviations from the original algorithm.

# Installation

The notebook can be viewed in github or using jupyter-notebook/-lab. To execute the code, some packages need to be installed to compute the eikonal equation .

Tested on Ubuntu 18.04. To run, switch to the repository directory and execute the following commands in a fresh python environment (e.g. virtual-env, or anaconda).

```bash
pip install -r requirements.txt
jupyter-notebook
```


# Limitations
While this repository is meant as a demonstration of the algorithm, several simplifications were made to make the code as short and concise as possible in contrast to the original paper:
- Projection directly from L-BFGS-B method instead of the Moreau-envelope
- A structured grid is used instead of a true FEM mesh
- \nabla \phi is estimated using a finite difference scheme
- Only the isotropic eikonal equation is considered, i.e. |\nabla \phi| = 1/v and solved using [skfmm](https://pythonhosted.org/scikit-fmm/)
- We assume quadrilateral basis functions inside the grid
- Contains neither the topological gradient, nor the ECG extension
- No handling of coalescing EASs

# Acknowledgements

If this work (or a part of it) helps you in your research, please consider acknowledging the github repository, or citing our paper.
The paper describing the method in detail is currently available on arXiv, but we will update the repository with the reference once it has
been published.

*TODO: Placeholder*
