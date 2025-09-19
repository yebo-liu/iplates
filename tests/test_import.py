def test_import_and_version():
import iplates


def test_euler_exports_exist():
    assert hasattr(iplates, 'compute_euler_pole')
    assert hasattr(iplates, 'compute_boundary_metrics')
    assert hasattr(iplates, 'BoundaryMetrics')

    assert hasattr(iplates, "__version__")
    assert isinstance(iplates.__version__, str)
    assert len(iplates.__version__) > 0


