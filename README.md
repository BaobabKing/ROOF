# ROOF
ROOF: Robustness-Oriented Optimal Factorization for Fluorescence Crosstalk Correction

This repository provides a MATLAB implementation of ROOF, a semi-blind crosstalk correction method for multi-band fluorescence microscopy. By incorporating a reference crosstalk matrix obtained from single-labeled control samples as a regularization prior, ROOF combines the strengths of experimental calibration with blind nonnegative matrix factorization, enabling robust, accurate, and efficient spectral unmixing.

---

The package contains four main subfolders:

- Algorithms/ – ROOF implementation, along with comparison methods, visualization tools, and performance evaluation utilities  
- Examples/ – Example scripts for running the algorithms  
- TestImages/ – Includes simulated and real fluorescence images  
- ExperimentResults/ – Sample outputs and benchmarking results  
