

from celest.coordinates.itrs import ITRS
from celest import units as u
from unittest import TestCase
import numpy as np
import os


class TestITRS(TestCase):

    def setUp(self):
        self.julian = np.random.rand(5)
        self.x = np.random.rand(5)
        self.y = np.random.rand(5)
        self.z = np.random.rand(5)
        self.unit = u.km
        self.valid_itrs = ITRS(self.julian, self.x, self.y, self.z, self.unit)

    def test_initialization(self):
        self.assertIsInstance(self.valid_itrs, ITRS)

    def test_parent_class_is_initialized(self):
        self.assertTrue(np.array_equal(self.julian, self.valid_itrs.time.data))
        self.assertEqual(u.jd2000, self.valid_itrs.time.unit)

        self.assertTrue(np.array_equal(self.x, self.valid_itrs.x.data))
        self.assertEqual(self.unit, self.valid_itrs.x.unit)

        self.assertTrue(np.array_equal(self.y, self.valid_itrs.y.data))
        self.assertEqual(self.unit, self.valid_itrs.y.unit)

        self.assertTrue(np.array_equal(self.z, self.valid_itrs.z.data))
        self.assertEqual(self.unit, self.valid_itrs.z.unit)

    def test_error_raised_with_different_length_inputs(self):
        julian = np.random.rand(1)
        x = np.random.rand(2)
        y = np.random.rand(3)
        z = np.random.rand(4)
        self.assertRaises(ValueError, ITRS, julian, x, y, z, self.unit)

    def test_error_raised_when_wrong_dimension_input_passed_in(self):
        julian = np.random.rand(1, 1)
        x = np.random.rand(2, 2)
        y = np.random.rand(3, 3)
        z = np.random.rand(4, 4)
        self.assertRaises(ValueError, ITRS, julian, x, y, z, self.unit)

    def test_save_to_text_file_saves_file(self):
        file_name = "itrs_test_file"
        self.valid_itrs.save_text_file(file_name)
        self.assertTrue(os.path.exists(file_name + ".txt"))
        if os.path.exists(file_name + ".txt"):
            os.remove(file_name + ".txt")
