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

    def __init__(
        self, experiment, model, fit_theta_bounds=(45.0, 110.0), model_kwargs={}
    ):
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
            `{}`.
        """

        if not isinstance(experiment, (str, Experiment)):
            raise ValueError(
                "Experiment should be either a string (filename) or Experiment class"
            )

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

        self.fit_theta_bounds = (np.min(fit_theta_bounds), np.max(fit_theta_bounds))
        self.exp_subset = (self.experiment.theta >= self.fit_theta_bounds[0]) & (
            self.experiment.theta <= self.fit_theta_bounds[1]
        )
        self.exp_theta_vals = self.experiment.theta[self.exp_subset]
        self.exp_int_vals = self.experiment.intensity[self.exp_subset]

        self.model_kwargs = dict(
            n_p=self.n_p,
            n_m=self.n_m,
            wavelength=self.wavelength,
            theta=self.experiment.theta,
        )

        if not "K_unit" in model_kwargs:
            model_kwargs["K_unit"] = experiment.K_unit

        if model_kwargs:
            self.model_kwargs = {**model_kwargs, **self.model_kwargs}
        # standard arguments that should be passed to each model instance

        self.fit_completed = False

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
        params = dict(zip(["prefactor", "diameter", "polydispersity"], fit_params))

        if self.debug:
            print(params)

        self.model = Model(
            d=params["diameter"], pd=params["polydispersity"], **self.model_kwargs
        )
        model_ints = (
            params["prefactor"]
            * self.model.intensity[np.isin(self.model.theta, self.exp_theta_vals)]
        )

        return 1.0 - (model_ints / self.exp_int_vals)

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

        self.fitter = least_squares(
            self._fit_function,
            start_params,
            bounds=[(0.0001, 200.0, 0.1), (1000.0, 2500.0, 20.0)],
            x_scale=(1.0e-3, 1.0, 1.0e-2),
            xtol=5e-4,
        )

        self.parameters = dict(
            zip(["prefactor", "diameter", "polydispersity"], self.fitter.x)
        )
        self.optimal_model = Model(
            d=self.parameters["diameter"],
            pd=self.parameters["polydispersity"],
            **self.model_kwargs,
        )

        self.fit_completed = True

        return self.optimal_model

    def save(self, filename=None, overwrite=True, extra_info={}, save_theta=False):
        """Saves a dictionary, to a pickle file, containing a `model_parameters`
        part which can directly be plugged into the `Model` class and a `fit_info`
        part which gives info on what file was fitted to, the theta range
        and the prefactor used for plotting model and experiment together.

        Parameters
        ----------
        filename : `str`, optional
            Filename to save the data as dictionary in pkl format, by default
            None, which results in saving it as [experiment_filename]_fit.pkl
        overwrite : `bool`, optional
            Whether or not to overwrite the file if it already exists, by
             default True
        extra_info : `dict`, optional
            Extra data in dictionary form to store in the `fit_info`, by default `{}`
        overwrite : `bool`, optional
            Whether or not to save the list of angles (theta), by default False

        Returns
        -------
        `bool`
            Whether or not it succeeded in saving the file
        """

        import os
        import pickle

        # set default filename if None is given
        if filename is None:
            filename = (
                str(self.experiment.filepath.resolve().absolute()).rsplit(".", 1)[0]
                + "_fit.pkl"
            )

        # save with correct extension
        if filename.rsplit(".", 1)[1] != "pkl":
            filename += ".pkl"

        # check if allowed to overwrite
        if not overwrite and os.path.exists(filename):
            raise FileExistsError(
                f"File {filename} already exists, set `overwrite` to `True` to"
                " overwrite this file or supply a different name"
            )

        if not self.fit_completed:
            raise AttributeError(
                "Data has not yet been fit, please call the `fit` function"
                " first before saving the results."
            )

        # construct dictionary for saving data
        # split into part that can directly be put into the Model class
        # and part with other information
        save_data = {}
        save_data["model_parameters"] = {**self.model_kwargs}
        if not save_theta:
            del save_data["model_parameters"]["theta"]
        save_data["model_parameters"]["d"] = self.parameters["diameter"]
        save_data["model_parameters"]["pd"] = self.parameters["polydispersity"]

        save_data["fit_info"] = {
            "filename": str(self.experiment.filepath.absolute()),
            "theta_range": self.fit_theta_bounds,
            "prefactor": self.parameters["prefactor"],
            **extra_info,
        }

        with open(filename, "wb") as datafile:
            pickle.dump(save_data, datafile)
            return True

        return False
