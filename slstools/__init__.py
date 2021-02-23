"""
Python package `slstools` contains several helper functions and classes
for [reading in data from the SCM SLS setup](#slstools.Experiment), [Miescat
 program](#slstools.Miescat) and for generating [models](#slstools.Model) 
for the scattering of spheres of given diameter and polydispersity.
The scattering calculations are based on [`miepython`](https://miepython.readthedocs.io/en/latest/).

Author: Roy Hoitink <L.D.Hoitink@uu.nl>
"""

from .experiment import Experiment
from .theory import Miescat
from .model import Model
from .fit import Fit

__all__ = ["Experiment", "Fit", "Model", "Miescat"]

# Do not build separate pages for submodules
__pdoc__ = {
    'experiment' : False,
    'fit': False,
    'theory' : False,
    'model' : False,
}