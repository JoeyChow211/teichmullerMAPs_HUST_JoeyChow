# Library Description

**Jacobian Determinant Validation for Lofted Hexahedral Macro-Elements** is a lightweight MATLAB code package for checking the local validity of three-dimensional lofted physical macro meshes. The library computes the center-point Jacobian determinant of each structured hexahedral macro-element and visualizes the determinant distribution using raw and median-normalized histograms.

The package was developed to support the quantitative validation of conformal-mapping-based lofted composite meshes. It is particularly useful for demonstrating that a generated physical macro mesh remains locally orientation-preserving and contains no element inversion. The default example evaluates a `16 x 10 x 20 = 3200` macro-element mesh and reports the minimum determinant, median determinant, normalized determinant range, and the number of non-positive elements.

This repository is intended for reproducible manuscript revision, peer-review verification, and academic reuse in mesh validation workflows for lofted or parameterized hexahedral domains.
