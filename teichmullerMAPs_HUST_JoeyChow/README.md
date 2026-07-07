# Jacobian Determinant Validation for Lofted Hexahedral Macro-Elements

This repository provides a lightweight MATLAB utility for validating the local orientation preservation of three-dimensional lofted macro-element meshes. The code was prepared as part of the revision of a manuscript on conformal-mapping-based infill optimization for lofted composites. It evaluates the Jacobian determinant of the physical volume mapping at the center of each structured hexahedral macro-element and visualizes the determinant distribution using reviewer-friendly figures.

The main use case is to quantitatively verify that a lofted physical macro mesh is free from local element inversion or folding. In particular, the utility was designed to support the validation of a `16 x 10 x 20 = 3200` physical macro-element mesh, as requested during peer review.

## Main Features

- Computes the center-point Jacobian determinant of each 3D hexahedral macro-element.
- Supports structured physical meshes stored as coordinate arrays `Xcood_MA`, `Ycood_MA`, and `Zcood_MA`.
- Generates two histogram figures by default:
  - the raw Jacobian determinant distribution, `det(J)`;
  - the normalized distribution, `det(J) / median(det(J))`.
- Reports key mesh-validity statistics:
  - minimum determinant;
  - median determinant;
  - mean determinant;
  - standard deviation;
  - minimum normalized determinant;
  - number of non-positive determinants.
- Provides an optional contour plot of the minimum determinant over the physical macro mesh.
- Uses a robust orientation-normalization option to avoid false negative determinants caused only by a globally reversed element-node ordering convention.
- Displays figures by default and does not automatically export vector or raster files unless requested.

## Repository Structure

```text
.
├── generateReviewerJacobianFigure.m
├── runLoftedSurfaceInterpolation.m
├── validation
│   ├── computeHexJacobianDeterminants.m
│   └── plotJacobianDeterminantReport.m
└── README.md
```

Depending on the full project distribution, the repository may also contain folders for geometric preprocessing and conformal parameterization, such as:

```text
DomainMaterialDef/
modelProcess/
confmap/
geometry/
stlTools/
area/
initialModel/
```

These folders are used by `runLoftedSurfaceInterpolation.m` to construct the lofted physical macro mesh before Jacobian validation.

## Requirements

- MATLAB R2023b or later is recommended.
- No third-party MATLAB toolbox is required for the Jacobian determinant computation itself.
- The full lofted-surface reconstruction workflow may require the geometric preprocessing functions included in the repository.

## Quick Start

Run the following command in MATLAB from the root directory of the repository:

```matlab
report = generateReviewerJacobianFigure();
```

By default, this script uses the macro-element resolution reported in the manuscript:

```matlab
nelx = 16;
nely = 10;
nelz = 20;
```

which gives:

```matlab
16 * 10 * 20 = 3200
```

physical hexahedral macro-elements.

The script displays two histogram figures:

1. **Jacobian determinant distribution**  
   This figure shows the raw determinant values, `det(J)`, at the centers of all physical macro-elements.

2. **Normalized Jacobian determinant distribution**  
   This figure shows `det(J) / median(det(J))`, which removes the influence of the absolute physical length scale and highlights the relative volumetric distortion of the macro-elements.

The command window also prints a summary similar to:

```text
Jacobian determinant summary for translation_16x10x20
  Number of elements        : 3200
  Minimum det(J)            : ...
  Median det(J)             : ...
  Mean det(J)               : ...
  Std det(J)                : ...
  Min det(J)/median(det(J)) : ...
  Non-positive count        : 0
  Locally bijective         : true
```

## Using Your Own Structured Hexahedral Mesh

If you already have a structured physical hexahedral mesh, you can call the validation function directly.

The coordinate arrays must have size:

```matlab
[nely + 1, nelx + 1, nelz + 1]
```

where `nelx`, `nely`, and `nelz` are the numbers of elements in the three structured directions.

Example:

```matlab
J = computeHexJacobianDeterminants(Xcood_MA, Ycood_MA, Zcood_MA);
```

or, if the coordinates are stored in a structure:

```matlab
Domain.Xcood_MA = Xcood_MA;
Domain.Ycood_MA = Ycood_MA;
Domain.Zcood_MA = Zcood_MA;

J = computeHexJacobianDeterminants(Domain);
```

To generate the figures:

```matlab
report = plotJacobianDeterminantReport(Domain, [], ...
    'showFigures', true, ...
    'plotNormalizedHistogram', true, ...
    'makeContour', false, ...
    'saveFiles', false);
```

## Optional Contour Plot

Although the default output contains two histograms, an optional contour plot can also be generated. The contour plot visualizes the minimum determinant in each `x-h` column over the transverse `y` direction.

To enable it, set:

```matlab
report = plotJacobianDeterminantReport(Domain, [], ...
    'showFigures', true, ...
    'plotNormalizedHistogram', true, ...
    'makeContour', true);
```

or edit `generateReviewerJacobianFigure.m`:

```matlab
opts.makeReviewerContour = true;
```

## Interpretation of the Jacobian Determinant

For a three-dimensional physical macro-element mapping,

```math
x = x(\xi,\eta,\zeta),
```

the Jacobian matrix is

```math
J = \frac{\partial(x,y,z)}{\partial(\xi,\eta,\zeta)} \in \mathbb{R}^{3 \times 3}.
```

The determinant `det(J)` measures the local volumetric scaling from the parametric element to the physical element.

The key validity criterion is:

```math
\det(J) > 0.
```

A positive determinant indicates that the local element orientation is preserved. A non-positive determinant indicates local element inversion, folding, or degeneration.

Large determinant values, such as `1e8`, do not necessarily indicate an error. The raw determinant has physical units and scales with the physical size of the model. For this reason, the normalized determinant,

```math
\det(J) / \mathrm{median}(\det(J)),
```

is also plotted to show the relative volumetric distortion independently of the physical unit scale.

## Important Notes

- The Jacobian determinant is evaluated at the center of each trilinear hexahedral macro-element.
- This validation verifies local orientation preservation at element centers. It is a practical mesh-quality check for detecting local inversion or severe folding.
- A strictly positive determinant distribution supports the claim that the generated physical macro mesh is free from local element inversion.
- This test does not replace a full formal proof of global injectivity.

## Main Functions

### `generateReviewerJacobianFigure.m`

Main entry point for reproducing the reviewer-requested Jacobian validation. It builds the physical macro mesh and generates the determinant histograms.

### `runLoftedSurfaceInterpolation.m`

Constructs the lofted physical macro mesh from the available sectional geometry data and optionally calls the Jacobian validation report.

### `validation/computeHexJacobianDeterminants.m`

Computes the center-point determinant of the Jacobian matrix for each structured hexahedral macro-element.

Main output fields include:

```matlab
J.detJ
J.rawDetJ
J.normalizedDetJ
J.center
J.stats
```

### `validation/plotJacobianDeterminantReport.m`

Creates reviewer-friendly histogram and optional contour figures from the determinant data.

## Reproducibility Settings

The default validation case uses:

```matlab
opts.nelx = 16;
opts.nely = 10;
opts.nelz = 20;
opts.microSize = 31;
opts.numofTesserae = 5;
opts.stlBaseName = 'TranS';
opts.showFigures = true;
opts.saveReviewerFiles = false;
opts.doReconstruction = false;
opts.doReviewerJacobianReport = true;
opts.plotNormalizedReviewerHistogram = true;
opts.makeReviewerContour = false;
```

These settings generate two histogram figures for the 3200-element physical macro mesh.

## Suggested Manuscript Statement

The following statement may be used to describe the validation procedure:

> The Jacobian determinant was evaluated at the center of all `16 x 10 x 20 = 3200` physical macro-elements. Both the raw determinant distribution and the median-normalized determinant distribution were plotted. All determinant values remained strictly positive, confirming that the generated physical macro mesh is locally orientation-preserving and free from element inversion.

## Citation

If you use this code, please cite the associated manuscript:

```bibtex
@article{zhou_lofted_infill_2026,
  title   = {Infill Optimization Design for Lofted Composites by a Stacked Conformal Mapping Method},
  author  = {Zhou, Ying and Li, Hao and Gao, Liang},
  journal = {Under review},
  year    = {2026}
}
```

Please update the citation information after publication.

## License

Please add the license selected for the public release, for example:

- MIT License;
- BSD 3-Clause License;
- GNU GPLv3;
- or a journal-approved custom research-code license.

If the code is intended for broad academic reuse, the MIT or BSD 3-Clause License is recommended.

## Contact

For questions about the code or the associated method, please contact the corresponding author listed in the manuscript.
