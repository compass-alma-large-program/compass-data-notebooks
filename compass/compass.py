"""
Utility functions for COMPASS data

"""

from astropy.io import fits
from astropy.table import QTable
from astropy import units as u


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


def maketable(
    source,
    datadir,
    settings=[1, 2, 3, 4, 5, 6, 7, 8, 9],
    spws=[25, 27, 29, 31],
    continuum=True,
):
    """
    Return a table with the RMS and synthesized beam...

    Parameters
    ----------
    source : str
        The name of the source.
    datadir : str
        The path to the data directory.
    settings : list, optional
        The settings to include in the table.
    spws : list, optional
        The spectral windows to include in the table.
    continuum : bool, optional
        If True, the continuum table is returned. If False, the line table is
        returned.

    Returns
    -------
    table : astropy.table.QTable

    """

    # TODO: Deal with missing values
    # TODO: Get data directory automatically

    t = QTable(
        names=(
            "Setting",
            "SPW",
            "Beam",
            "RMS",
            "Thermal Noise",
            "Number of pixels",
            "Pixel size",
            "Image size",
        ),
        dtype=[int, int, str, float, float, int, float, float],
        units=[None, None, None, "mJy / beam", "mJy / beam", None, "marcsec", "arcsec"],
    )

    for setting in settings:
        for spw in spws:
            # Read the FITS file.
            if continuum:
                filename = (
                    f"{datadir}/{source}/cubes/{source}-set{setting}-spw{spw}-cont.fits"
                )
            else:
                filename = f"{datadir}/{source}/cubes/{source}-set{setting}-spw{spw}-lines.fits"
            try:
                hdul = fits.open(filename)
            except FileNotFoundError:
                break  # next setting

            # Get the beam
            bmaj = hdul[0].header["BMAJ"] * u.deg
            bmin = hdul[0].header["BMIN"] * u.deg
            bpa = hdul[0].header["BPA"] * u.deg
            bmaj = bmaj.to(u.arcsec).value
            bmin = bmin.to(u.arcsec).value
            bpa = bpa.to(u.deg).value
            beam_str = f"{bmaj:3.2f}'' x {bmin:3.2f}'' ({bpa:3.0f}°)"

            # Get the noise
            rms = hdul[0].header["NOISEMEA"] * u.Jy / u.beam
            thnoise = hdul[0].header["NOISETHE"] * u.Jy / u.beam  # thermal noise

            # Get number of pixels and the pixel size
            npix = hdul[0].header["NAXIS1"]
            pixsize = abs(hdul[0].header["CDELT1"] * u.deg)
            imsize = npix * pixsize

            # Add values to the table
            t.add_row(
                (
                    setting,
                    spw,
                    beam_str,
                    rms,
                    thnoise,
                    npix,
                    pixsize,
                    imsize,
                )
            )

            hdul.close()

    # Format columns
    t["Pixel size"].info.format = ".0f"
    t["Image size"].info.format = ".0f"
    if continuum:
        t["RMS"] = t["RMS"].to(u.uJy / u.beam)
        t["Thermal Noise"] = t["Thermal Noise"].to(u.uJy / u.beam)
        t["RMS"].info.format = ".0f"
        t["Thermal Noise"].info.format = ".0f"
    else:
        t["Thermal Noise"].info.format = ".2f"
        t["RMS"].info.format = ".2f"

    return t
