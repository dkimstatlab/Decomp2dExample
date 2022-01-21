# Decomp2dExample
Decomposition examples for two-dimensional image

The subjects of image processing include feature extraction, compression, denoising, image enhancement, and restoration. 
Many software have been developed for these studies. However, few studies have focused on decomposition in the literature, 
even though image decomposition is essential for image processing. We provide an R solution termed "Decomp2d", 
tailored to the viewpoint of image decomposition based on existing R packages.

The "Decomp2d" solution can be installed through the following commands.
```
> install.packages("devtools")
> devtools::install_github("dkimstatlab/Decomp2d")
> library(Decomp2d)
```
Or the package can be installed using package archive file "Decomp2d_0.6.0.tar.gz".

We illustrate the usage of the "Decomp2d" solution and describe the procedure using an example image and D11 texture.

- Decomp2d-manual.pdf : manual of R package "Decomp2d" 
- Decomp_D11.R : R code for decomposition of D11 texture 
- Decomp_examples.R : R code for decomposition of synthethic two-dimensional image 
