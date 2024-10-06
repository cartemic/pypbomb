import numpy as np

from pypbomb.data import bolt
from pypbomb.data.bolt import ExternalThreadAnsi, InternalThreadAnsi, ThreadType


def test_available_thread_sizes():
    assert "2-6" in bolt.available_thread_sizes(ThreadType.external, "2A")
    assert "2-6" not in bolt.available_thread_sizes(ThreadType.external, "3A")
    assert "2-6" in bolt.available_thread_sizes(ThreadType.internal, "2B")
    assert "2-6" not in bolt.available_thread_sizes(ThreadType.internal, "3B")


def test_available_thread_classes():
    external_2a = bolt.available_thread_sizes(ThreadType.external, "2A")
    external_3a = bolt.available_thread_sizes(ThreadType.external, "3A")
    assert "0-80" in external_2a
    assert "0-80" in external_3a
    assert "2-6" in external_2a
    assert "2-6" not in external_3a

    internal_2b = bolt.available_thread_sizes(ThreadType.internal, "2B")
    internal_3b = bolt.available_thread_sizes(ThreadType.internal, "3B")
    assert "0-80" in internal_2b
    assert "0-80" in internal_3b
    assert "2-6" in internal_2b
    assert "2-6" not in internal_3b


def test_internal_thread_properties():
    result = bolt.dimensions(ThreadType.internal, thread_size="9/16-20", thread_class="2B")
    assert isinstance(result, InternalThreadAnsi)
    assert result.thread_size == "9/16-20"
    assert result.thread_class == "2B"
    assert result.threads_per_inch == 20
    assert np.isclose(result.diameter.basic, 0.0142875)  # 0.5625 in -> m
    assert np.isclose(result.diameter.pitch.min, 0.013462)  # 0.53 in -> m
    assert np.isclose(result.diameter.pitch.max, 0.0136017)  # 0.5355 in -> m
    assert np.isclose(result.diameter.minor.min, 0.0129032)  # 0.508 in -> m
    assert np.isclose(result.diameter.minor.max, 0.013208)  # 0.52 in -> m
    assert np.isclose(result.diameter.major, 0.0142875)  # 0.5625 in -> m


def test_external_thread_properties():
    result = bolt.dimensions(ThreadType.external, thread_size="9/16-20", thread_class="2A")
    assert isinstance(result, ExternalThreadAnsi)
    assert result.thread_size == "9/16-20"
    assert result.thread_class == "2A"
    assert result.threads_per_inch == 20
    assert np.isclose(result.diameter.basic, 0.0142875)  # 0.5625 in -> m
    assert np.isclose(result.diameter.pitch.min, 0.0133223)  # 0.5245 in -> m
    assert np.isclose(result.diameter.pitch.max, 0.01342898)  # 0.5287 in -> m
    assert np.isclose(result.diameter.major.min, 0.01425448)  # 0.5612 in -> m
    assert np.isclose(result.diameter.major.max, 0.01404874)  # 0.5531 in -> m
    assert np.isclose(result.diameter.minor, 0.01274318)  # 0.5017 in -> m
