"""
Utility functions for COMPASS data

"""

from astropy.io import fits


def extract_spectrum_from_cube(cube, position):
    """
    Extract a spectrum from a cube.

    Parameters
    ----------
    cube : SpectralCube
        The cube to extract from.
    position : SkyCoord
        The position to extract from.

    Returns
    -------
    spectra : SpectralCube
        The extracted spectra.

    """

    # Get the pixel position.
    dec, ra = map(int, position.to_pixel(cube.wcs))

    # Extract the spectrum.
    try:
        spectrum = cube[:, dec, ra]
    except IndexError:
        raise IndexError("Position is outside the cube.")

    return spectrum


def fixradesys(filename):
    """
    Fix the RADESYS keyword in the header of a FITS file.

    Parameters
    ----------
    filename : str
        The name of the FITS file to fix.

    """

    with fits.open(filename, mode="update") as hdul:
        hdr = hdul[0].header
        if hdr["RADESYS"] == "FK5":
            hdr["RADESYS"] = "ICRS"
            hdul.flush()
            hdul.close()

    return
