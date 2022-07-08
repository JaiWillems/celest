

from celest.file_save import TextFileWriter, _process_file_name
from unittest import TestCase, skip
import numpy as np
import os


class TestTextFileWriter(TestCase):

    def setUp(self):
        self.file_name = "file_name"
        self.processed_file_name = _process_file_name(self.file_name)
        self.header = "header"
        self.writer = TextFileWriter(self.file_name, self.header)

        self.subheading = "subheading"
        self.parameters = [
            ["label 1", "value 1"],
            ["label 2", "value 2"]
        ]
        self.data = [
            ["data 1", np.random.rand(5)],
            ["data 2", np.random.rand(5)]
        ]
        self.data_format = "%.5f"

    def tearDown(self):
        if os.path.exists(self.processed_file_name):
            os.remove(self.processed_file_name)

    def test_initialize(self):
        self.assertIsInstance(self.writer, TextFileWriter)

    def test_process_file_name_ends_in_txt_when_not_in_file_name(self):
        file_name = "test_file"
        self.assertEqual(file_name + ".txt", _process_file_name(file_name))

    def test_process_file_name_ends_in_txt_when_in_file_name(self):
        file_name = "test_file.txt"
        self.assertEqual(file_name, _process_file_name(file_name))

    def test_adding_layer_adds_layer(self):
        self.assertEqual(0, len(self.writer._layers))
        self._add_test_layer()
        self.assertEqual(1, len(self.writer._layers))
        self._add_test_layer()
        self.assertEqual(2, len(self.writer._layers))

    def _add_test_layer(self):
        self.writer.add_layer(self.subheading, self.parameters, self.data,
                              self.data_format)

    def test_save_creates_file(self):
        self._add_test_layer()
        self._add_test_layer()
        self.writer.save()
        self.assertTrue(os.path.exists(self.processed_file_name))
