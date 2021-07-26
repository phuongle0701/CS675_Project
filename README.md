# CS675_Project: A survey of Krylov Subspace methods for nonsymmetric linear systems.
Project in CS675 - Computational Linear Algebra

## Direction: 

1. The folder 'Krylov_subspace_method' contains implementation of `GMRES`, `BiCG`, and `IDR(s)` methods. 

2. The folder 'PoissonPDE' builds the finite difference scheme for {1D,2D,3D} for Poisson PDE equations with 2 heat sources temperature. The `Laplacian` returns the matrix A corresponds to the lhs of finite difference scheme and matrix b of source terms. 

3. The folder 'Image_Denoising' builds the PDE-constrained equation to denoise the camera-man image (tif). 




## References:

```
@book{saad_book,
author = {Saad, Y.},
title = {Iterative Methods for Sparse Linear Systems},
year = {2003},
isbn = {0898715342},
publisher = {Society for Industrial and Applied Mathematics},
address = {USA},
edition = {2nd}
}


@article{Saad1986GMRESAG,
  title={GMRES: a generalized minimal residual algorithm for solving nonsymmetric linear systems},
  author={Y. Saad and M. Schultz},
  journal={Siam Journal on Scientific and Statistical Computing},
  year={1986},
  volume={7},
  pages={856-869}
}

@article{idr_method,
author = {Sonneveld, Peter and Gijzen, M.B.},
year = {2008},
month = {01},
pages = {1035-1062},
title = {IDR( s ): A Family of Simple and Fast Algorithms for Solving Large Nonsymmetric Systems of Linear Equations},
volume = {31},
journal = {SIAM J. Scientific Computing},
doi = {10.1137/070685804}
}
```

