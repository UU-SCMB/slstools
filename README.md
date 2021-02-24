# SLS python helpers

Opiniated set of helper functions to handle measurments from SCM's self-built static light scattering (SLS) setup and fit models to it.

## Info

Author: Roy Hoitink <L.D.Hoitink@uu.nl>

## Installation

Installation of this package can be done via pip. 
To install the newest version of this package, run the following command:

```bash
pip install git+https://github.com/rhoitink/slstools --upgrade
```

This ensures that the package and its dependencies are installed.

## Documentation

The documentation can be found at [rhoitink.github.io/slstools](https://rhoitink.github.io/slstools).

## Example usage

Several helper classes are available for loading experiments (`Experiment` class), creating models (`Model` class) and fitting a model to an experiment (`Fit` class). 
Below an example is given on how to fit a model to an experiment.

```python
# import the helpers and libaries
from slstools import Experiment, Model, Fit
import matplotlib.pyplot as plt

# load the experiment that is saved in the file 'Sample01.sls'
experiment = Experiment("Sample01.sls", K_unit="nm")

# Create an initial model of 1000nm (diameter) spheres with a polydispersity of 5%
# Refractive index medium: 1.333 (water) and particle: 1.4345 (n-hexadecane)
# Make sure the given diameter rougly matches the expected diameter, as this will improve the fitting (and its speed)
model = Model(d=1000, pd=5, n_p=1.4345, n_m=1.333)

# Load both model and experiment into the Fit class
# Fitting will occur for scattering angles between 45 and 110 degrees
fit = Fit(experiment, model, fit_theta_bounds=(45.0, 110.0),model_kwargs=dict(K_unit=experiment.K_unit))

# Obtain the optimal model after fitting
optimal_model = fit.fit()

# plot the experimental data and the fit
# scale the optimal_model with the found prefactor such that the lines overlap
plt.plot(optimal_model.K, fit.parameters["prefactor"]*optimal_model.intensity, label=f"Fit: d={optimal_model.diameter:.0f}nm ({optimal_model.polydispersity:.0f}%)")

# plot experimental data as comparison
plt.plot(experiment.K, experiment.intensity, 'k.', alpha=.3, label="Experimental data")

# setting labels, scales and legend
plt.xlabel(r"K (%s$^{-1}$)" % optimal_model.K_unit)
plt.ylabel("Normalised intensity (a.u.)")
plt.yscale("log")
plt.legend()
plt.show()
```
