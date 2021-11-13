"""Testing module for the `Interpolation` class."""


from celest.core.interpolation import _interpolate
from unittest import TestCase
import numpy as np
import unittest


class TestInterpolation(TestCase):

    def test_interp(self):
        """Test `Interpolation._interp`."""

        data = np.linspace(0, 99, 100)

        for factor in [0, 1, 2]:

            # factor = {0, 1, 2}, dt=0, indices=None
            out_data = np.linspace(0, 99, factor * 100)
            calc_out = _interpolate(data=data, factor=factor)

            for i in range(calc_out.size):
                with self.subTest(i=i):
                    self.assertAlmostEqual(out_data[i], calc_out[i], delta=0.001)

            for dt in [0, 1, 2]:

                if factor == 1:
                    alpha = 1
                else:
                    alpha = 0

                # factor = {0, 1, 2}, dt={0, 1, 2}, indices=[[10, 11, ..., 20]]
                indices = np.array([[i for i in range(10, 21)]])

                out_1 = np.linspace(0, 9 - dt, 10 - dt)
                num = factor * (10 + 2 * dt + alpha)
                out_2 = np.linspace(10 - dt, 20 + dt, num)
                out_3 = np.linspace(21 + dt, 99, 79 - dt)
                out_data = np.concatenate((out_1, out_2, out_3))
                calc_out = _interpolate(data=data, factor=factor, dt=dt, indices=indices)

                for i in range(calc_out.size):
                    with self.subTest(i=i):
                        self.assertAlmostEqual(out_data[i], calc_out[i], delta=0.001)

                # factor = {0, 1, 2}, dt={0, 1, 2}, indices=[[10, 11, ..., 20],
                # [40, 41, ..., 50]]
                indices = np.array([[i for i in range(10, 21)], [i for i in range(40, 51)]])

                out_1 = np.linspace(0, 9 - dt, 10 - dt)
                num = factor * (10 + 2 * dt + alpha)
                out_2 = np.linspace(10 - dt, 20 + dt, num)
                out_3 = np.linspace(21 + dt, 39 - dt, 19 - 2 * dt)
                num = factor * (10 + 2 * dt + alpha)
                out_4 = np.linspace(40 - dt, 50 + dt, num)
                out_5 = np.linspace(51 + dt, 99, 49 - dt)
                out_data = np.concatenate((out_1, out_2, out_3, out_4, out_5))
                calc_out = _interpolate(data=data, factor=factor, dt=dt, indices=indices)

                for i in range(calc_out.size):
                    with self.subTest(i=i):
                        self.assertAlmostEqual(out_data[i], calc_out[i], delta=0.001)


if __name__ == "__main__":
    unittest.main()
