# test_profile.py

"""Tests for profile.py"""

from conftest import *
from src import profile
from unittest import mock

# Test retrieval and encoding of recordings from sequence and config file
def test_get_codes_from_file(mock_fastq, mock_settings, mock_codes_from_fastq):
    with mock.patch("builtins.open", side_effect=mock.mock_open(read_data=mock_fastq)) as mock_input_f, mock.patch("yaml.safe_load", return_value=mock_settings) as mock_config_f:
        recordings = profile.get_recordings("mock.fastq", "mock.yaml")
        assert recordings == mock_codes_from_fastq

# Test retrieval and encoding of recordings from sequence file with invalid recordings
def test_get_codes_from_file_with_invalid_recordings(mock_fastq_with_invalid_recordings, mock_settings, mock_codes_from_fastq_with_invalid_recordings):
    with mock.patch("builtins.open", side_effect=mock.mock_open(read_data=mock_fastq_with_invalid_recordings)) as mock_input_f, mock.patch("yaml.safe_load", return_value=mock_settings) as mock_config_f:
        recordings = profile.get_recordings("mock.fastq", "mock.yaml")
        assert recordings == mock_codes_from_fastq_with_invalid_recordings