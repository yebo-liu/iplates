def test_import_and_version():
    import iplates

    assert hasattr(iplates, "__version__")
    assert isinstance(iplates.__version__, str)
    assert len(iplates.__version__) > 0


