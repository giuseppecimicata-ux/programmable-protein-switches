# programmable-protein-switches
Python tools for quantitative analysis and thermodynamic modeling of programmable protein conformational switches
# Programmable Protein Switches

Python tools for quantitative analysis and thermodynamic modeling of programmable protein conformational switches.

## Overview

This repository contains Python analysis and modeling workflows developed to characterize programmable protein conformational switches controlled by competitive protein–protein interactions.

The experimental system uses a binder to modulate the conformational state of a protein switch and a competing switch-release molecule to sequester the binder and reverse the conformational transition.

The computational workflow integrates experimental titrations, nonlinear parameter estimation, mechanistic equilibrium modeling, two-dimensional switching landscapes, and quantitative comparison between model predictions and experimental measurements.

## Computational workflow

The analysis proceeds from independently measured experimental parameters to prediction and validation of the two-dimensional switching landscape:

1. **Spectral characterization** — visualization of fluorescence titration measurements and experimental variability.

2. **Binder titration** — nonlinear fitting of normalized titration data to estimate an apparent dissociation constant and Hill coefficient.

3. **Switch-release countertitration** — mechanistic modeling of competitive binder sequestration to estimate the apparent interaction strength of the switch-release molecule.

4. **Reversible switching cycles** — visualization of repeated binder/switch-release cycles to assess reversible molecular actuation.

5. **Experimental 2D switching landscape** — visualization of normalized switching activity across combinations of binder and switch-release concentrations.

6. **Thermodynamic prediction** — prediction of the two-dimensional switching landscape from independently measured interaction parameters using equilibrium models of increasing mechanistic detail.

7. **Model validation** — quantitative comparison of experimentally measured and model-predicted switching landscapes.

## Repository structure

```text

scripts/
    plot_titration_spectra.py
    fit_binder_titration.py
    fit_switch_release_countertitration.py
    plot_reversible_switching_cycles.py
    plot_experimental_2d_landscape.py
    predict_2d_landscape_simple_model.py
    predict_2d_landscape_coupled_model.py
    validate_2d_switching_model.py

docs/
    model_description.md
```

## Scripts

### `plot_titration_spectra.py`

Plots fluorescence spectra from protein-switch titration experiments. Error bars represent the standard deviation of triplicate measurements.

### `fit_binder_titration.py`

Fits normalized binder-titration data using a Hill model to estimate an apparent dissociation constant and Hill coefficient. Experimental variability is shown as the standard deviation of triplicate measurements.

### `fit_switch_release_countertitration.py`

Fits switch-release countertitration data using a mechanistic 1:1 binding model that explicitly accounts for binder sequestration and the resulting change in free binder concentration.

### `plot_reversible_switching_cycles.py`

Plots FRET measurements across repeated binder/switch-release cycles to visualize reversible protein-switch actuation. Error bars represent the standard deviation of triplicate measurements.

### `plot_experimental_2d_landscape.py`

Generates a two-dimensional response landscape by mapping normalized switching activity across combinations of binder and switch-release concentrations.

### `predict_2d_landscape_simple_model.py`

Predicts the two-dimensional switching landscape from independently measured apparent interaction constants using a simplified equilibrium model.

### `predict_2d_landscape_coupled_model.py`

Implements a coupled equilibrium model of switch, binder, and switch-release interactions to generate a mechanistically refined prediction of the two-dimensional switching landscape.

### `validate_2d_switching_model.py`

Compares experimentally measured and model-predicted switching landscapes using a parity plot and reports the coefficient of determination (R²) as a measure of agreement.

## Experimental variability

Where indicated, error bars represent the standard deviation of triplicate experimental measurements.

## Model assumptions

The analysis uses apparent interaction constants derived from independent titration experiments. Apparent EC50 values are interpreted as approximate apparent dissociation constants under the binding assumptions described in the associated manuscript.

The thermodynamic models make simplifying assumptions appropriate to the experimental system. Detailed equations and model assumptions are documented in `docs/model_description.md`.

## Requirements

The analysis uses Python and standard scientific-computing libraries including:

* NumPy
* pandas
* SciPy
* Matplotlib
* openpyxl

Exact dependencies and versions are provided in `requirements.txt`.

## Associated manuscript

**[Manuscript title]**

[Authors]

[Journal / preprint / DOI when available]

## Author

**Giuseppe Cimicata, PhD**
University of California, San Francisco

Research areas: protein engineering, structural biology, quantitative biophysics, and computational protein design.
