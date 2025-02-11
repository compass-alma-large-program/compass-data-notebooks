"""
Functions to apply the primary beam correction

"""

import numpy as np
from scipy.special import j1
from astropy import units as u
from astropy import constants as const
from matplotlib import pyplot as plt


def alma_primary_beam(theta, freq):
    """
    Compute the primary beam of the ALMA 12m telescope.

    This function assumes that the primary beam is that of uniformly
    illuminated 10.7 m dish with 0.75 m blockage.

    Parameters
    ----------
    theta : `float` or array-like
        The angular distance from the beam center, in arcsec.
    freq : `float`
        The frequency in GHz.

    Returns
    -------
    beam : `float` or array-like
        The primary beam.

    Notes
    -----
    See: https://help.almascience.org/kb/articles/how-do-i-model-the-alma-primary-beam-and-how-can-i-use-that-model-to-obtain-the-sensitivity-pr

    """

    D = 10.7 * u.m  # diameter of the dish
    d = 0.75 * u.m  # diameter of the subreflector

    theta *= u.arcsec
    freq *= u.GHz
    wl = const.c / freq  # wavelength

    D = D.to(u.m).value
    d = d.to(u.m).value
    theta = theta.to(u.rad).value
    freq = freq.to(u.Hz).value
    wl = wl.to(u.m).value

    # Intensity pattern of a uniformely illuminated circular aperture
    # with a diameter D, with a central obscuration of diameter d.
    alpha = (d / D) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        beam = (
            2 * j1(np.pi * D * theta / wl) / (np.pi * D * theta / wl)
            - alpha * 2 * j1(np.pi * d * theta / wl) / (np.pi * d * theta / wl)
        ) ** 2 / (1 - alpha) ** 2
    beam = np.nan_to_num(beam, nan=1.0)  # Replace NaN at theta=0 with 1

    return beam


def alma_primary_beam_gaussian(theta, freq):
    """
    Compute a Gaussian a Gaussian approximation of the ALMA beam.

    This function assumes that the primary beam is a Gaussian with a
    FWHM of 1.13 λ/D, with D = 12 m.

    Parameters
    ----------
    theta : `float` or array-like
        The angular distance from the beam center, in arcsec.
    freq : `float`
        The frequency in GHz.

    Returns: array-like
        The primary beam
    """

    D = 12 * u.m  # Antenna diameter

    theta *= u.arcsec
    freq *= u.GHz
    wl = const.c / freq  # wavelength

    D = D.to(u.m).value
    theta = theta.to(u.rad).value
    freq = freq.to(u.Hz).value
    wl = wl.to(u.m).value
    fwhm = 1.13 * wl / D
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to sigma

    return np.exp(-(theta**2) / (2 * sigma**2))


def apply_primary_beam_correction(cube, threshold=0.2):
    """
    Apply the primary beam correction to a cube

    Parameters
    ----------
     cube : `SpectralCube`
        The spectral cube to correct.
     threshold : `float`
        The threshold to apply to the primary beam mask (fraction of
        the peak).

    """

    # Compute the primary beam
    freq = (cube.wcs.wcs.restfrq * u.Hz).to(u.GHz).value
    velo, dec, ra = cube.world[0, :, :]
    ra0, dec0 = cube.header["RA"] * u.deg, cube.header["DEC"] * u.deg
    theta = (
        np.sqrt(((ra - ra0) * np.cos(dec0)) ** 2 + (dec - dec0) ** 2).to(u.arcsec).value
    )  # offset from the phase center
    pb = alma_primary_beam(theta, freq)

    # Apply the primary beam correction
    cube = cube.with_mask(pb > threshold)
    # FIXME: This would be more efficient if we could use this, but it does not work
    # See: https://github.com/radio-astro-tools/spectral-cube/issues/939
    # cube.apply_numpy_function(np.divide, cube, pb, reduce=False)
    cube.allow_huge_operations = True
    cube = cube / pb

    return cube


def plot_primary_beam():
    """
    Plot the ALMA primary beam at 300 GHz.
    """

    freq = 300  # GHz
    theta = np.linspace(-60, 60, 1000)

    plt.plot(
        theta, alma_primary_beam(theta, freq), label="10.7-m dish with 0.75m blockage"
    )
    plt.plot(
        theta, alma_primary_beam_gaussian(theta, freq), label="Gaussian approximation"
    )
    plt.xlabel(r"$\theta$ (arcsec)")
    plt.ylabel("Primary Beam")
    plt.legend()
    plt.show()
