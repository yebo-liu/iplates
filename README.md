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



### Euler pole calculation

This package estimates an Euler pole that moves an oceanic plate toward the dominant subduction directions, with angular speed proportional to the subduction fraction of the plate boundary.

Notation:
- Unit Cartesian for a lat/lon point:
```
x(phi, lambda) = (cos(phi)cos(lambda), cos(phi)sin(lambda), sin(phi))
```
- Great-circle midpoint for unit vectors `a, b`:
```
m_hat = normalize(a + b)   # skip if a ≈ -b (antipodal)
```
- Segment weight `w_i` is proportional to segment length.

1) Boundary centroid (approximate):
```
c = normalize( sum_i w_i * m_hat_i )
```

2) Subduction-driven surface direction at `c`:
- For each subduction segment, average its midpoints to get unit `s`.
- Plane normal: `n = normalize( c × s )`.
- Tangent direction at `c` toward subduction: `t = normalize( n × c )`.
- Weighted mean direction:
```
T = sum_j w_j * t_j
t_hat = normalize( T )
```

3) Euler pole axis direction (so that `(omega × c) ≈ t_hat`):
```
omega_hat = normalize( c × t_hat )
```

4) Angular speed and total rotation:
Let `f_sub ∈ [0,1]` be the subduction fraction (subduction length / total boundary length). Given user coefficient `k` (deg/Myr) and age span `Δt` (Myr):
```
omega_deg_per_ma     = f_sub * k
total_rotation_angle = omega_deg_per_ma * Δt
```

API example:
```python
from iplates import compute_euler_pole

pole = compute_euler_pole(model, plate_id=701, age_ma=0.0, age_span=10.0, velocity_coefficient=1.5)
# {
#   'euler_lat': ..., 'euler_lon': ...,
#   'omega_deg_per_ma': f_subduction * velocity_coefficient,
#   'total_rotation_angle_deg': omega_deg_per_ma * age_span,
#   'subduction_fraction': f_subduction, 'metrics': {...}
# }
```

