## Introduction
IBFest is the research code for ion beam figuring (IBF) system for synchrotron X-ray mirrors designed and developed by the Optical Metrology and Fabrication Group at National Synchrotron Light Source II (NSLS-II), NY, US. 

## Implemented algorithms
- Rectangular surface error map simulation using 2D Lengendre polynomials
- 2D Beam Removal Function (BRF) fittintg and learning
- Dwell time calculation algorithms
  - Fourier domain methods
    - [x] Bayesian iterative method [4]
    - [x] Fourier transform + Inverse filtering [2]
    - [x] Robust Iterative Fourier Trasform-based dwell time Algorithm (RIFTA) [8, 9]
  - Matrix-based methods
    - [x] Truncated SVD (TSVD) [3]
    - [x] LSQR [1, 5]
    - [x] Constrained Linear Least Squares (CLLS) + Coarse-to-Fine scheme [6, 7]
- Thresholded inverse filtering assisted by Nelder-Mead Simplex algorithm [9]
- High-performance 2D convolution using FFT

## Usage
**Note:**
  - The CLLS algorithm applies 'active-set' algorithm to solve the CLLS equations. This function has been removed from MATLAB since 2016b. This algorithm is preferred since it converges faster than the other two. 
  - All the units used in the code are metres unless otherwise specified.
  - To properly run the TSVD, make sure the computer has at least 16GB RAM since SVD consumes a lot memory.

The common arguments that are required for IBFest dwell time calculation include:
  - Beam Removal Function (BRF): the BRF can come from either the measurement or model, by choosing ```avg``` or ```model```, respectively. If 'model' is chosen, the parameters for a 2D Gaussian should be set, includeing the Peak Removal Rate (PRR) ```A```, the ```Sigma```, the diameter ```d```, and the centers ```u```. If 'avg' is chosen, ```X_brf```, ```Y_brf```, and ```Z_brf``` should be provided. 
  - Z_to_remove: the desired height to be removed. 
  - ca_range: the range of the Clear Aperture (CA), which is a struct contains ```x_s```, ```y_s```, ```x_e```, and ```y_e```, which are the start and end coordinates (units are meters) of the CA. 
  - dw_range: the range of the DWell grid (DW), which should be larger than ca_range at least the radius of the BRF on each side.

## Reference

[1] [Carnal, C. L., Egert, C. M., & Hylton, K. W. (1992, December). Advanced matrix-based algorithm for ion-beam milling of optical components. In Current Developments in Optical Design and Optical Engineering II (Vol. 1752, pp. 54-62). International Society for Optics and Photonics.](https://doi.org/10.1117/12.130719)

[2] [Wilson, S. R., & McNeil, J. R. (1987, January). Neutral ion beam figuring of large optical surfaces. In Current Developments in Optical Engineering II (Vol. 818, pp. 320-324). International Society for Optics and Photonics.](https://doi.org/10.1117/12.978903)

[3] [Zhou, L., Dai, Y. F., Xie, X. H., Jiao, C. J., & Li, S. Y. (2007). Model and method to determine dwell time in ion beam figuring. Nanotechnol. Precis. Eng., 5(8–9), 107-112.](http://en.cnki.com.cn/Article_en/CJFDTotal-NMJM200702009.htm)

[4] [Jiao, C., Li, S., & Xie, X. (2009). Algorithm for ion beam figuring of low-gradient mirrors. Applied Optics, 48(21), 4090-4096.](https://doi.org/10.1364/AO.48.004090)

[5] [Wu, J. F., Lu, Z. W., Zhang, H. X., & Wang, T. S. (2009). Dwell time algorithm in ion beam figuring. Applied optics, 48(20), 3930-3937.](https://doi.org/10.1364/AO.48.003930)

[6] [Wang, T., Huang, L., Vescovi, M., Kuhne, D., Tayabaly, K., Bouet, N., & Idir, M. (2019). Study on an effective one-dimensional ion-beam figuring method. Optics express, 27(11), 15368-15381.](https://doi.org/10.1364/OE.27.015368)

[7] [Wang, T., Huang, L., Vescovi, M., Kuhne, D., Tayabaly, K., Bouet, N., & Idir, M. (2019, September). One-dimensional ion-beam figuring solution from Brookhaven National Laboratory. In Advances in Metrology for X-Ray and EUV Optics VIII (Vol. 11109, p. 1110909). International Society for Optics and Photonics.](https://doi.org/10.1117/12.2526074)

[8] [Wang, T., Huang, L., Tayabaly, K., & Idir, M. (2019, November). Study on the performances of dwell time algorithms in ion beam figuring. In Optifab 2019 (Vol. 11175, p. 111750M). International Society for Optics and Photonics.](https://doi.org/10.1117/12.2536869)

[9] T. Wang, L. Huang, H. Kang, H. Choi, D. W. Kim, K. Tayabaly, andM. Idir, “Rifta: a robust iterative fourier transform-based dwell time algorithm for ion beam figuring,” Sci. Reports, Accepted (2020)

[10] [T. Wang, L. Huang, Y. Zhu, M. Vescovi, D. Khune, H. Kang, H. Choi, D. Kim, K. Tayabaly, N. Bouet, and M. Idir, "Development of a position–velocity–time-modulated two-dimensional ion beam figuring system for synchrotron x-ray mirror fabrication," Appl. Opt.  59, 3306-3314 (2020).](https://doi.org/10.1364/AO.389010)
