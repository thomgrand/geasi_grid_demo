# GEASI Demo
This


# Installation

Prerequisites
-------------
Tested on Ubuntu 18.04. To run, switch to the repository directory and execute the following commands in a fresh python environment (e.g. virtual-env, or anaconda).

```bash
pip install -r requirements.txt
jupyter-notebook
```


Limitations
----------
While this repository is meant as a demonstration of the algorithm, several simplifications were made to make the code as short and concise as possible:
- Projection directly from L-BFGS-B method instead of the Moreau-envelope
- A structured grid is used instead of a true FEM mesh
- $\nabla \phi$ is estimated using a finite difference scheme
- Only the isotropic eikonal equation is considered, i.e. $\Vert \nabla \phi \Vert = \frac{1}{v}$
- We assume quadrilateral basis functions inside the grid
- Contains neither the topological gradient, nor the ECG extension
- No handling of coalescing EASs

# Acknowledgements

If this works helps you in your research, please consider acknowledging the github repository, or citing our paper.
The paper describing the method in detail is currently available on arXiv, but we will update the repository with the reference once it has
been published.

*TODO: Placeholder*
