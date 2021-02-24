from .experiment import Experiment
from .model import Model
from scipy.optimize import least_squares
import numpy as np

class Fit:
    """Class for performing fitting a model to experimental data.

    ## Example usage
    ```python
    from slstools import Experiment, Model, Fit
    experiment = Experiment("Sample01.sls", K_unit="nm")
    model = Model(d=1000, pd=5, n_p=1.4345, n_m=1.333)

    fit = Fit(experiment, model, model_kwargs=dict(K_unit="nm"))

    optimal_model = fit.fit()
    ```
    """

    def __init__(self, experiment, model, fit_theta_bounds=(45.,110.), model_kwargs=None):
        """Initialize fitting a given `Model` to given `Experiment`

        Parameters
        ----------
        experiment : `str` or `Experiment`
            The experiment that is used for fitting, if a `str` is given
            this will be passed as filename to the `Experiment` class.
        model : `Model`
            Intial model that will be used for the fit, make sure to 
            choose a sensible starting diameter and polydispersity as it helps
            to speed up the fitting. Also use the correct refractive indices
            for the particle and medium.
        fit_theta_bounds : `tuple`, optional
            Range of `theta` (scattering angle) that will be used for fitting
            data, by default `(45.,110.)`.
        model_kwargs : `dict`, optional
            Extra keyword arguments to pass onto the final model, by default
            `None`.
        """

        if type(experiment) not in [str, Experiment]:
            raise ValueError("Experiment should be either a string (filename) or Experiment class")

        if type(experiment) == str:
            self.experiment = Experiment(experiment)
        else:
            self.experiment = experiment
        
        self.model = model

        self.n_p = self.model.n_particle
        self.n_m = self.model.n_medium
        self.d0 = self.model.diameter
        self.pd0 = self.model.polydispersity
        self.wavelength = self.model.lambda0

        self.fit_theta_bounds = (np.min(fit_theta_bounds),np.max(fit_theta_bounds))
        self.exp_subset = (self.experiment.theta >= self.fit_theta_bounds[0]) & (self.experiment.theta <= self.fit_theta_bounds[1])
        self.exp_theta_vals = self.experiment.theta[self.exp_subset]
        self.exp_int_vals = self.experiment.intensity[self.exp_subset]

        self.model_kwargs = dict(n_p = self.n_p, n_m=self.n_m, wavelength=self.wavelength, theta=self.experiment.theta)
        if model_kwargs:
            self.model_kwargs = {**model_kwargs, **self.model_kwargs}
        # standard arguments that should be passed to each model instance
        

    def _fit_function(self, fit_params):
        """Internal function that is used as cost function for fitting

        Parameters
        ----------
        fit_params : `list`
            List of fitting parameters in the form of 
            `[prefactor, diameter, polydispersity]`.

        Returns
        -------
        np.ndarray
            List of ratios of model_intensity/experiment_intensity as measure
            of the quality of the fit.
        """
        params = dict(zip(['prefactor','diameter','polydispersity'], fit_params))
        
        if self.debug:
            print(params)
        
        self.model = Model(d=params['diameter'], pd=params['polydispersity'], **self.model_kwargs)
        model_ints = params['prefactor']*self.model.intensity[np.isin(self.model.theta, self.exp_theta_vals)]
        
        return 1.0-(model_ints/self.exp_int_vals)

    def fit(self, debug=False):
        """Function to perform the actual fitting of the model to the experiment

        Parameters
        ----------
        debug : `bool`, optional
            Whether or not to output intermediate fitting parameters,
            by default `False`.

        Returns
        -------
        `Model`
            The Model with parameters that best fit the experiment
            Note: the intensities should be multiplied with 
            `self.optimal_params[0]` to be accurate to the experiment's
            intensities.
        """
        self.debug = debug

        start_params = [1.0, self.d0, self.pd0]

        self.fitter = least_squares(self._fit_function, start_params, bounds=[(.0001,200.,.1),(1000.,2500.,20.)],\
            x_scale=(1.0e-3, 1.0, 1.0e-2), xtol=1e-3)

        self.parameters = dict(zip(["prefactor", "diameter", "polydispersity"],self.fitter.x))
        self.optimal_model = Model(d=self.parameters["diameter"], pd=self.parameters["polydispersity"], **self.model_kwargs)

        return self.optimal_model