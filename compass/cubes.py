"""
Tools to manipulate line datacubes

"""

import matplotlib.axes as maxes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.visualization.wcsaxes import add_beam


def extract_line(cube, restfreq, vmin, vmax):
    """
    Extract a line from a datacube

    Parameters
    ----------
    cube : SpectralCube
        The datacube
    restfreq : Quantity
        The rest frequency of the line
    vmin : Quantity
        The minimum velocity to extract
    vmax : Quantity
        The maximum velocity to extract

    Returns
    -------
    line : SpectralCube
        The extracted line
    """

    return cube.with_spectral_unit(
        u.km / u.s, rest_value=restfreq, velocity_convention="radio"
    ).spectral_slab(vmin, vmax)


def noise_K(cube):
    """
    Return the noise in K

    Parameters
    ----------
    cube : SpectralCube

    Returns
    -------
    noise : float
        The noise in K
    """

    return (cube.header["NOISEMEA"] * u.Jy / u.beam).to(
        u.K, equivalencies=cube.beam.jtok_equiv(cube.header["RESTFRQ"])
    )


def zeroth_order_moment(
    cube,
    fig,
    map_center,
    map_size,
    vmin=None,
    vmax=None,
    levels=None,
    subplot=111,
    **kwargs,
):
    """
    Compute the 0th moment and map it

    This function returns a WCSAxes instance containing a map of the
    zeroth-order moment of a datacube.

    Parameters
    ----------
    cube : SpectralCube
        Name of the FITS file
    fig : `matplotlib.figure.Figure`
        Figure instance in which the map will be draw
    map_center : `astropy.coordinates.SkyCoord`
        Center of the map
    map_size : Quantity
        Size of the map
    vmin : float, optional
        Minimum value to display
    vmax : float, optional
        Maximum value to display
    levels : array-like, optional
        Contour levels to overplot on the map
    subplot : int, optional
        3-digit integer that give the subplot position

    Returns
    -------
    ax : `astropy.visualization.wcsaxes.WCSAxes`
        Axes instance containing the map_size

    cb : `matplotlib.colorbar.Colorbar`
        Colorbar instance corresponding to the map

    """

    # Compute the 0th order moment
    moment0 = cube.moment(order=0).to(u.K * u.km / u.s)

    # Make a color map with 0th order moment
    ax = fig.add_subplot(subplot, projection=moment0.wcs, **kwargs)
    im = ax.imshow(moment0.data, vmin=0, vmax=moment0.max().value * 1.1)

    # Overplot the 0th order moment contours
    if levels is not None:
        ax.contour(moment0.data, levels=levels, colors="white", linewidths=1)

    # Set limits
    map_center_pix = map_center.to_pixel(wcs=moment0.wcs)
    pix_scale = proj_plane_pixel_scales(ax.wcs)
    degree_per_pixel = (pix_scale[0] * pix_scale[1]) ** 0.5 * u.degree
    ax.set_xlim(
        map_center_pix[0] - 0.5 * map_size / degree_per_pixel,
        map_center_pix[0] + 0.5 * map_size / degree_per_pixel,
    )
    ax.set_ylim(
        map_center_pix[1] - 0.5 * map_size / degree_per_pixel,
        map_center_pix[1] + 0.5 * map_size / degree_per_pixel,
    )

    # Set the range of the color bar
    if vmin is not None:
        im.set_clim(vmin=vmin)
    if vmax is not None:
        im.set_clim(vmax=vmax)

    # Add a color bar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="8%", pad=0, axes_class=maxes.Axes)
    cb = fig.colorbar(im, cax=cax)

    # Show the beam size
    add_beam(
        ax,
        moment0.header,
        corner="bottom right",
        frame=True,
        borderpad=0,
        pad=1,
        color="black",
    )

    return ax, cb
