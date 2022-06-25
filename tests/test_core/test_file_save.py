

import os
from celest.core.file_save import _save_data_as_txt, get_parameter_string
from unittest import TestCase


class TestSaveDataAsTxt(TestCase):

    def setUp(self):

        self.file_name = "test_file"

    def tearDown(self):

        if os.path.exists(self.file_name + ".txt"):
            os.remove(self.file_name + ".txt")

    def test_save_data_creates_file(self):
        _save_data_as_txt(self.file_name)
        self.assertTrue(os.path.exists(self.file_name + ".txt"))

    def test_get_parameter_string_for_non_final_parameter(self):
        parameter_string = get_parameter_string(0, 2, "label", "value")
        self.assertEqual("label: value\n", parameter_string)

    def test_get_parameter_string_for_final_parameter(self):
        parameter_string = get_parameter_string(0, 1, "label", "value")
        self.assertEqual("label: value", parameter_string)
