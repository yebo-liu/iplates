### iplates

Personal utilities for plate reconstruction and plotting. This is a lightweight, private package intended for my own workflows.

### Install (editable for development)

From this directory (`Jupyter/iplates`):

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .[dev]
```

### Run tests

```bash
pytest -q
```

### Usage

After installing, import as usual:

```python
import iplates
print(iplates.__version__)

# Optionally access modules (these may require external deps):
# from iplates.reconstruction import resolve_topos, plot_boundaries
```

### Notes

- This project is for personal use; API may change without notice.
- Heavy GIS dependencies (e.g., pygplates, geopandas) are not pinned here. Install them in your environment as needed.


