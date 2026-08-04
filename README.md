# Doppler-Based EKF Correction of SGP4/TLE Predictions

## Overview

This repository contains MATLAB code developed for a masters project on improving low Earth orbit (LEO) satellite predictions. The project investigates whether Doppler shift measurements from ground stations can be used to correct trajectory predictions produced from Two-Line Element (TLE) data using the Simplified General Perturbations 4 (SGP4) model.

SGP4/TLE propagation is computationally efficient, but its prediction accuracy generally decreases as the propagation time from the TLE epoch increases. The implemented method first estimates and, when necessary, removes a significant temporal offset between the predicted and reference observations. An extended Kalman filter (EKF) then uses Doppler residuals to update selected orbital parameters.

## Method

The main processing sequence is:

1. Read the TLE data and propagate the satellite orbit using SGP4.
2. Calculate or import ground-station observations, including Doppler shift, azimuth, and elevation.
3. Estimate the temporal offset between the predicted and reference observations.
4. Remove the temporal offset when its magnitude is greater than one minute.
5. Apply the EKF to correct selected orbital parameters using Doppler measurements.
6. Evaluate the correction using Doppler, azimuth, elevation, error plots, and root mean square error (RMSE).

The EKF state vector used in the project is based on selected TLE/SGP4 parameters:

$$
\mathbf{x} = [M,\ \Omega,\ \omega,\ i,\ e,\ n,\ B^*]^T,
$$

where $M$ is mean anomaly, $\Omega$ is the right ascension of the ascending node, $\omega$ is the argument of perigee, $i$ is inclination, $e$ is eccentricity, $n$ is mean motion, and $B^*$ is the drag term.

## Main Features

- SGP4 propagation from TLE data
- Ground-station Doppler, azimuth, and elevation calculations
- Satellite visibility filtering based on elevation
- Temporal-offset estimation and alignment
- Doppler-based EKF correction of orbital parameters
- Support for single- and multiple-ground-station experiments
- Error analysis and RMSE calculation
- Optional additive white Gaussian noise (AWGN) for simulated measurements

## Requirements

- MATLAB
- The MATLAB functions and data files included in this repository
- MATLAB Aerospace Toolbox may be required by scripts that use toolbox-specific satellite or coordinate functions

## Usage

1. Clone or download this repository.
2. Open the project folder in MATLAB.
3. Add the project and its subfolders to the MATLAB path.
4. Check the TLE file paths, observation-data file paths, satellite selection, ground-station coordinates, and EKF parameters in the script you intend to run.
5. Run the relevant experiment script.

Input files may include TLE files and reference or simulated observation data. The scripts generate predicted and corrected observations together with comparison plots and error metrics.

## Original SGP4 Code and Modifications

This project includes code derived from the following open-source MATLAB implementation:

> Meysam Mahooti (2026). [SGP4](https://www.mathworks.com/matlabcentral/fileexchange/62013-sgp4), MATLAB Central File Exchange. Retrieved August 4, 2026.

The original SGP4 submission calculates orbital state vectors for near-Earth satellites from TLE data. It is distributed under the BSD 3-Clause License.

The modifications and extensions made for this thesis project include:

- integration of ground-station Doppler observations;
- estimation and removal of temporal offsets;
- implementation and tuning of an extended Kalman filter;
- correction of selected orbital parameters;
- single- and multiple-ground-station test configurations; and
- additional plotting, visibility filtering, and error analysis.

The original author retains copyright in the original SGP4 code. Copyright in the modifications and newly added code belongs to Jingxun Li.

## License

This repository is distributed under the BSD 3-Clause License. The original copyright notice for Meysam Mahooti is retained, together with the copyright notice for the modifications. See [`LICENSE`](LICENSE) for details.

## References

1. Hoots, F. R., and Roehrich, R. L. (1980). *Models for Propagation of NORAD Element Sets*. Spacetrack Report No. 3, U.S. Air Force Aerospace Defense Command.
2. Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications* (4th ed.). Microcosm Press.
3. Mahooti, M. (2026). [SGP4](https://www.mathworks.com/matlabcentral/fileexchange/62013-sgp4), MATLAB Central File Exchange.

