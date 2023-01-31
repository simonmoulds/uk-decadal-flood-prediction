#!/usr/bin/env python3

from decimal import Decimal

import iris
import numpy as np

# Default fill-value.
_MDI = 1e+20

# Stock cube - global grid extents (degrees).
_LAT_MIN = -90.0
_LAT_MAX = 90.0
_LAT_RANGE = _LAT_MAX - _LAT_MIN
_LON_MIN = 0.0
_LON_MAX = 360.0
_LON_RANGE = _LON_MAX - _LON_MIN

def _generate_cube_from_dimcoords(latdata, londata, circular: bool = False):
    """Generate cube from lat/lon points.
    Parameters
    ----------
    latdata : np.ndarray
        List of latitudes.
    londata : np.ndarray
        List of longitudes.
    circular : bool
        Wrap longitudes around the full great circle. Bounds will not be
        generated for circular coordinates.
    Returns
    -------
    :class:`~iris.cube.Cube`
    """
    lats = iris.coords.DimCoord(latdata,
                                standard_name='latitude',
                                units='degrees_north',
                                var_name='lat',
                                circular=circular)

    lons = iris.coords.DimCoord(londata,
                                standard_name='longitude',
                                units='degrees_east',
                                var_name='lon',
                                circular=circular)

    if not circular:
        # cannot guess bounds for wrapped coordinates
        lats.guess_bounds()
        lons.guess_bounds()

    # Construct the resultant stock cube, with dummy data.
    shape = (latdata.size, londata.size)
    dummy = np.empty(shape, dtype=np.dtype('int8'))
    coords_spec = [(lats, 0), (lons, 1)]
    cube = iris.cube.Cube(dummy, dim_coords_and_dims=coords_spec)

    return cube


def _spec_to_latlonvals(*, start_latitude: float, end_latitude: float,
                        step_latitude: float, start_longitude: float,
                        end_longitude: float, step_longitude: float) -> tuple:
    """Define lat/lon values from spec.
    Create a regional cube starting defined by the target specification.
    The latitude must be between -90 and +90. The longitude is not bounded, but
    wraps around the full great circle.
    Parameters
    ----------
    start_latitude : float
        Latitude value of the first grid cell center (start point). The grid
        includes this value.
    end_latitude : float
        Latitude value of the last grid cell center (end point). The grid
        includes this value only if it falls on a grid point. Otherwise, it
        cuts off at the previous value.
    step_latitude : float
        Latitude distance between the centers of two neighbouring cells.
    start_longitude : float
        Latitude value of the first grid cell center (start point). The grid
        includes this value.
    end_longitude : float
        Longitude value of the last grid cell center (end point). The grid
        includes this value only if it falls on a grid point. Otherwise, it
        cuts off at the previous value.
    step_longitude : float
        Longitude distance between the centers of two neighbouring cells.
    Returns
    -------
    xvals : np.array
        List of longitudes
    yvals : np.array
        List of latitudes
    """
    if step_latitude == 0:
        raise ValueError('Latitude step cannot be 0, '
                         f'got step_latitude={step_latitude}.')

    if step_longitude == 0:
        raise ValueError('Longitude step cannot be 0, '
                         f'got step_longitude={step_longitude}.')

    if (start_latitude < _LAT_MIN) or (end_latitude > _LAT_MAX):
        raise ValueError(
            f'Latitude values must lie between {_LAT_MIN}:{_LAT_MAX}, '
            f'got start_latitude={start_latitude}:end_latitude={end_latitude}.'
        )

    def get_points(start, stop, step):
        """Calculate grid points."""
        # use Decimal to avoid floating point errors
        num = int(Decimal(stop - start) // Decimal(str(step)))
        stop = start + num * step
        return np.linspace(start, stop, num + 1)

    latitudes = get_points(start_latitude, end_latitude, step_latitude)
    longitudes = get_points(start_longitude, end_longitude, step_longitude)

    return latitudes, longitudes


def regional_stock_cube(spec: dict):
    """Create a regional stock cube.
    Returns
    -------
    :class:`~iris.cube.Cube`.
    """
    latdata, londata = _spec_to_latlonvals(**spec)

    cube = _generate_cube_from_dimcoords(latdata=latdata,
                                         londata=londata,
                                         circular=True)

    def add_bounds_from_step(coord, step):
        """Calculate bounds from the given step."""
        bound = step / 2
        points = coord.points
        coord.bounds = np.vstack((points - bound, points + bound)).T

    add_bounds_from_step(cube.coord('latitude'), spec['step_latitude'])
    add_bounds_from_step(cube.coord('longitude'), spec['step_longitude'])

    return cube
