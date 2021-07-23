"""
Python package `slstools` contains several helper functions and classes
for [reading in data from the SCM SLS setup](#slstools.Experiment), [Miescat
 program](#slstools.Miescat), generating [models](#slstools.Model) 
for the scattering of spheres of given diameter and polydispersity and
[fitting](#slstools.Fit) those models to the experimental data.
The scattering calculations are based on [`miepython`](https://miepython.readthedocs.io/en/latest/).

Author: Roy Hoitink <L.D.Hoitink@uu.nl>
"""
from .experiment import Experiment
from .theory import Miescat
from .model import Model
from .fit import Fit

__version__ = "0.2.0"

__all__ = ["Experiment", "Fit", "Model", "Miescat"]

# Do not build separate pages for submodules
__pdoc__ = {
    "experiment": False,
    "fit": False,
    "theory": False,
    "model": False,
}
