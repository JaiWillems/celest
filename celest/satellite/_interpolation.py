"""General interpolation functionality."""


from celest.core.decorators import set_module
from scipy.interpolate import interp1d
import numpy as np


@set_module('celest.satellite')
class Interpolation(object):
    """General interpolation object.

    The `Interpolation` class allows for general interpolation without an
    independent data array.

    Methods
    -------
    _interp(data, factor, dt, indices)
        Interpolate `data` using the given parameters.
    """

    def _interp(self, data: np.array, factor: int=0, dt: int=0, indices:
                 np.array=None) -> np.array:
        """Interpolate `data` using the given parameters.

        Parameters
        ----------
        data : np.array
            Array of shape (i, j) containing columns of similar data.
        factor : int
            The factor increase in the number of sample points.
        dt : int
            The number of points laying beyond the `indices` range to
            interpolate within.
        indices : np.array
            Array of array containing the interpolation indices.

        Returns
        -------
        np.array
            The `data` parameter interpolated with the provided parameters.

        Note
        ----
        The dt parameter is only relevent if an indices parameter is given.

        Examples
        --------
        >>> interp = Interpolation()
        >>> data = np.array([[1, 6], [2, 7], [3, 8], [4, 9], [5, 10]])
        >>> interp(data=data, factor=2, dt=1, indices=np.array([[2, 3]]))
        array([[ 1.        ,  6.        ],
               [ 2.        ,  7.        ],
               [ 2.42857143,  7.42857143],
               [ 2.85714286,  7.85714286],
               [ 3.28571429,  8.28571429],
               [ 3.71428571,  8.71428571],
               [ 4.14285714,  9.14285714],
               [ 4.57142857,  9.57142857],
               [ 5.        , 10.        ]])
        """

        j = data.shape[0]

        indep = np.linspace(0, j - 1, j)
        interp = interp1d(x=indep, y=data, axis=0, kind="cubic")

        if type(indices) == type(None):

            indep_new = np.linspace(0, j - 1, j * factor)

        else:

            indep_new = np.linspace(0, j - 1, j)

            for reg in np.flipud(indices):

                extension = np.arange(reg[-1] + 1, reg[-1] + 1 + dt, 1)
                reg = np.append(reg, extension)
                reg = np.insert(reg, 0, np.arange(reg[0] - dt, reg[0], 1))

                min_I, max_I = reg[0], reg[-1]

                if min_I != max_I:
                    num_ind = int(factor * (max_I - min_I + 1))
                    reg_indep_new = np.linspace(min_I, max_I, num_ind)
                else:
                    reg_indep_new = reg

                delete_ind = reg[np.where(reg <= j - 1)[0]]
                indep_new = np.delete(indep_new, delete_ind)

                insert_vals = reg_indep_new[np.where(reg_indep_new <= j - 1)[0]]
                indep_new = np.insert(indep_new, min_I, insert_vals)

        data_new = interp(indep_new)

        return data_new
