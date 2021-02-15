import numpy as np
import miepython as mp


class Model:
    """Helper class to create a model for scattering spheres of given diameter
    (and polydispersity).

    Attributes
    ----------
    K : np.ndarray
        Numpy array containing the scattering angles in units of `K_unit`.
    intensity : np.ndarray
        Numpy array containing scattering intensities for each angle
        (perpendicular polarisation).
    """

    def __init__(
        self,
        n_p,
        d=1000.0,
        n_m=1.333,
        wavelength=632.8,
        theta=np.arange(0.0, 181.0, 0.5),
        pd=0.0,
        K_unit="m",
        normalise=True,
    ):
        """Initialisation of the class, assigns data to class attributes and
        calculates model from given parameters.

        Parameters
        ----------
        n_p : float|complex
            (complex) refractive index of the particle.
        d : float, optional
            Diameter of the particle in `nm`, by default `1000.0`.
        n_m : float, optional
            Refractive index of the medium, by default `1.333` (water).
        wavelength : float, optional
            Vacuum wavelength of the laser in `nm`,
            by default 632.8 (HeNe laser).
        theta : np.ndarray, optional
            List of angles to calculate intensities for,
            by default `np.arange(0.0,181.0,.5)`
        pd : float, optional
            Polydispersity of the system in % (e.g. `5` for 5% p.d.),
            by default `0.0`.
        K_unit : str, optional
            (inverse) unit of scattering angle K, by default `m`.
        normalise : bool, optional
            Whether or not to normalise the intensities to `1.0`,
            by default `True`.
        """
        if K_unit not in ["m", "nm", "um", "µm"]:
            raise ValueError(
                f"The unit of K should be m, nm or µm. The given value '{K_unit}' is not allowed."
            )

        self.n_particle = n_p
        self.n_medium = n_m
        self.n_ratio = self.n_particle / self.n_medium
        self.lambda0 = wavelength
        self.lambdam = wavelength / self.n_medium

        self.theta = theta
        if type(self.theta) is not np.ndarray:
            self.theta = np.array(self.theta)

        self.diameter = d

        self.K = (
            4
            * np.pi
            * self.n_medium
            / (self.lambda0 * 1e-9)
            * np.sin(self.theta / 2.0 * np.pi / 180)
        )
        self.K_unit = K_unit

        if self.K_unit == "nm":
            self.K /= 1e9
        elif self.K_unit == "um" or self.K_unit == "µm":
            self.K_unit = "µm"
            self.K /= 1e6

        if pd > 0.0:
            # created gaussian distribution with given polydispersity

            self.intensity = np.zeros(self.theta.shape)
            sigma = self.diameter * (pd / 100.0)
            dist_width = 3 * sigma
            dist_diams = np.linspace(
                self.diameter - dist_width, self.diameter + dist_width, 9
            )
            dist_weights = (
                1
                / (sigma * np.sqrt(2 * np.pi))
                * np.exp(-((dist_diams - self.diameter) ** 2.0) / (2.0 * sigma ** 2.0))
            )

            for weight, diam in zip(dist_weights, dist_diams):
                self.intensity += weight * self._single_scatterer(diam)
        else:
            self.intensity = self._single_scatterer()

        if normalise:
            self.intensity /= self.intensity.max()

    def _single_scatterer(self, diameter=None):
        """Function to generate scattering intensities for given diameter

        Parameters
        ----------
        diameter : float, optional
            Diameter of the sphere in `nm`. If `None`, the class attribute
            `diameter` will be used, by default `None`

        Returns
        -------
        np.ndarray
            Numpy array containing the scattering intensities (normalised for
            scattering volume, perpendicular polarisation).
        """
        if not diameter:
            diameter = self.diameter
        mu = np.cos(self.theta * np.pi / 180.0)
        size_param = 2 * np.pi / self.lambdam * diameter / 2
        geometric_cross_section = np.pi * (diameter * 1e-9 / 2.0) ** 2

        qext, qsca, *_ = mp.ez_mie(
            self.n_particle, diameter, self.lambda0, self.n_medium
        )
        ipar, iper = mp.ez_intensities(
            self.n_particle, diameter, self.lambda0, mu, self.n_medium
        )

        sigma_sca = geometric_cross_section * qsca * iper

        return sigma_sca
