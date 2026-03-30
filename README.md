# LIM and CSLIM MATLAB Implementation

This repository contains MATLAB code for constructing:

- Linear Inverse Models (LIM)
- Cyclostationary Linear Inverse Models (CSLIM)

from oceanic variability data.

The implementation uses principal component (PC) time series derived from sea surface temperature (SST) and sea surface height (SSH) to estimate linear dynamical operators from lagged covariance statistics. Both time-invariant (LIM) and seasonally varying (CSLIM) models are included, along with example forecasts.

---

## Repository Structure

```
├── code/
   ├── LIM_CSLIM_construction.m
   └── nearestSPD.m

├── data/
   ├── EOFs.nc
   ├── LIM_CSLIM_operators.nc
   └── state_vector_data.nc

```
---

## Methods Overview

The code:

- Constructs a reduced-order state vector from SST and SSH principal components  
- Computes covariance and lagged covariance matrices  
- Estimates linear dynamical operators governing system evolution  

- Builds a **Linear Inverse Model (LIM)** using a single, time-invariant operator  

- Builds a **Cyclostationary LIM (CSLIM)** using month-dependent operators that vary seasonally  

- Computes associated noise covariance matrices for each model  
- Performs example forecasts using the learned operators  

---

## References

The repository includes the `nearestSPD` algorithm (D’Errico, 2026), which ensures covariance matrices are symmetric positive definite. This is required for numerical stability and for use with multivariate normal sampling methods (e.g., `mvnrnd` in MATLAB).

D’Errico, J. (2026). nearestSPD. MATLAB Central File Exchange.  
https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

---

## Data

All required input data are included.

The NetCDF files provide:
- Principal component time series
- Empirical orthogonal functions (EOFs)
- Precomputed operator matrices

These allow the model construction and example forecasts to be reproduced directly.

---

## Usage

1. Open MATLAB  
2. Navigate to the repository directory  
3. Run:

LIM_CSLIM_construction

Make sure the `code/` and `data/` folders are accessible in your MATLAB path.

---

## Requirements

- MATLAB  
- NetCDF support

---

## Software Archive

A snapshot of this repository is archived on Zenodo:

https://doi.org/10.5281/zenodo.19258346

---

## Authors  

Daniel Vimont  
University of Wisconsin–Madison  

Alexandra Vizcarra  
University of Wisconsin–Madison  

Jack Zweifel  
University of Wisconsin–Madison

---
