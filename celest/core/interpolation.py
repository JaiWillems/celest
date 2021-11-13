"""General interpolation functionality.

The module contains the functions to allow for general interpolation of column
data without an independent data array.
"""


from celest.core.decorators import set_module
from scipy.interpolate import interp1d
import numpy as np


def _interpolate(data: np.ndarray, factor: int=1, dt: int=0, indices:
                np.ndarray=None) -> np.ndarray:
    """Interpolate `data` per the given parameters.

    This method takes a multicolumn data array and interpolates it along
    each column, in the same manner, using the interpolation parameters.

    Parameters
    ----------
    data : np.ndarray
        Array of shape (i, j) containing columns of similar information. Note
        that the data must be two dimensional.
    factor : int, optional
        The factor increase in the number of sample points of interpolated
        regions.
    dt : int, optional
        A positive integer defining the number of points laying beyond the
        `indices` range to interpolate within.
    indices : np.ndarray, optional
        Array of non-empty arrays containing the interpolation indices in a
        monotonically increasing fashion.

    Returns
    -------
    np.ndarray
        The `data` parameter interpolated per the given parameters.

    Notes
    -----
    The `_interp` method was designed to allow for greater precision data
    in regions where it is important. As a result, it allows for the
    interpolation of full data sets or only where may be necessary.

    The `dt` parameter is only relavent if an `indices` parameter is given.

    Examples
    --------
    >>> data_1 = np.linspace(0, 9, 10).reshape((-1, 1))
    >>> data_2 = np.linspace(5, 14, 10).reshape((-1, 1)))
    >>> data = np.concatenate((data_1, data_2, axis=1)
    >>> ind = np.array([[2, 4]]
    >>> Interpolation()._interp(data=data, factor=2, dt=1, indices=ind))
    np.array([[ 0.          5.        ]
              [ 1.          6.        ]
              [ 1.57142857  6.57142857]
              [ 2.14285714  7.14285714]
              [ 2.71428571  7.71428571]
              [ 3.28571429  8.28571429]
              [ 3.85714286  8.85714286]
              [ 4.42857143  9.42857143]
              [ 5.         10.        ]
              [ 6.         11.        ]
              [ 7.         12.        ]
              [ 8.         13.        ]
              [ 9.         14.        ]])
    """

    j = data.shape[0]

    indep = np.linspace(0, j - 1, j, dtype=int)
    interp = interp1d(x=indep, y=data, axis=0, kind="cubic")

    if (indices is None) or (indices.size == 0):

        indep_new = np.linspace(0, j - 1, j * factor)

    else:

        indep_new = np.copy(indep)

        for reg in np.flipud(indices):

            min_I, max_I = int(reg[0] - dt), int(reg[-1] + dt)

            if min_I < 0:
                min_I = 0

            if max_I > indep[-1]:
                max_I = indep[-1]

            alpha = 1 if factor == 1 else 0

            left_reg = indep_new[:min_I]

            new_reg = np.linspace(min_I, max_I, factor * (max_I - min_I + alpha))
            right_reg = indep_new[max_I + 1:]

            indep_new = np.concatenate((left_reg, new_reg, right_reg))

    data_new = interp(indep_new)

    return data_new
