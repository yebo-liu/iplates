import math
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

try:
    import pygplates
except ImportError:
    pygplates = None


EARTH_RADIUS_KM = 6371.0


@dataclass
class BoundaryMetrics:
    plate_id: int
    age_ma: float
    total_length_km: float
    subduction_length_km: float
    ridge_length_km: float
    transform_length_km: float
    other_length_km: float
    subduction_fraction: float
    centroid_lat: Optional[float]
    centroid_lon: Optional[float]


def _ensure_pygplates():
    if pygplates is None:
        raise ImportError("pygplates is required for iplates.euler but is not installed in this environment.")


def _polyline_length_km(geometry) -> float:
    points = getattr(geometry, 'get_points', lambda: None)()
    if not points:
        return 0.0
    length_rad = 0.0
    for i in range(1, len(points)):
        arc = pygplates.GreatCircleArc(points[i - 1], points[i])
        length_rad += arc.get_arc_length()
    return length_rad * EARTH_RADIUS_KM


def _polyline_midpoints(points) -> List[pygplates.PointOnSphere]:
    mids: List[pygplates.PointOnSphere] = []
    if not points or len(points) < 2:
        return mids
    for i in range(1, len(points)):
        a = points[i - 1]
        b = points[i]
        a_xyz = _point_to_xyz(a)
        b_xyz = _point_to_xyz(b)
        m_xyz = _normalize(_add(a_xyz, b_xyz))
        if m_xyz == (0.0, 0.0, 0.0):
            # Opposite points; skip midpoint
            continue
        m_lat, m_lon = _xyz_to_lat_lon(m_xyz)
        mids.append(pygplates.PointOnSphere(m_lat, m_lon))
    return mids


def _point_to_xyz(point: pygplates.PointOnSphere) -> Tuple[float, float, float]:
    lat_deg, lon_deg = point.to_lat_lon()
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    x = math.cos(lat) * math.cos(lon)
    y = math.cos(lat) * math.sin(lon)
    z = math.sin(lat)
    return (x, y, z)


def _xyz_to_lat_lon(xyz: Tuple[float, float, float]) -> Tuple[float, float]:
    x, y, z = xyz
    r = math.sqrt(x * x + y * y + z * z)
    if r == 0:
        return (None, None)
    x, y, z = x / r, y / r, z / r
    lat = math.degrees(math.asin(z))
    lon = math.degrees(math.atan2(y, x))
    return (lat, lon)


def _normalize(v: Tuple[float, float, float]) -> Tuple[float, float, float]:
    x, y, z = v
    n = math.sqrt(x * x + y * y + z * z)
    if n == 0:
        return (0.0, 0.0, 0.0)
    return (x / n, y / n, z / n)


def _cross(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    ax, ay, az = a
    bx, by, bz = b
    return (ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx)


def _add(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    ax, ay, az = a
    bx, by, bz = b
    return (ax + bx, ay + by, az + bz)


def _scale(a: Tuple[float, float, float], s: float) -> Tuple[float, float, float]:
    ax, ay, az = a
    return (ax * s, ay * s, az * s)


def _dot(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    ax, ay, az = a
    bx, by, bz = b
    return ax * bx + ay * by + az * bz


def _rotate_about_axis(v: Tuple[float, float, float],
                       axis: Tuple[float, float, float],
                       angle_rad: float) -> Tuple[float, float, float]:
    """
    Rodrigues' rotation formula to rotate unit vector v about unit axis by angle_rad.
    """
    vx, vy, vz = v
    ax, ay, az = _normalize(axis)
    cos_t = math.cos(angle_rad)
    sin_t = math.sin(angle_rad)
    # v_parallel = (a·v) a
    a_dot_v = ax * vx + ay * vy + az * vz
    v_parallel = (ax * a_dot_v, ay * a_dot_v, az * a_dot_v)
    # v_perp = v - v_parallel
    v_perp = (vx - v_parallel[0], vy - v_parallel[1], vz - v_parallel[2])
    # a × v
    a_cross_v = _cross((ax, ay, az), (vx, vy, vz))
    rx = v_parallel[0] + v_perp[0] * cos_t + a_cross_v[0] * sin_t
    ry = v_parallel[1] + v_perp[1] * cos_t + a_cross_v[1] * sin_t
    rz = v_parallel[2] + v_perp[2] * cos_t + a_cross_v[2] * sin_t
    return _normalize((rx, ry, rz))


def _resolve_topology_for_plate(topology_model, age_ma: float, plate_id: int):
    snapshot = topology_model.topological_snapshot(age_ma)
    for topo in snapshot.get_resolved_topologies():
        try:
            pid = topo.get_resolved_feature().get_reconstruction_plate_id()
        except Exception:
            try:
                pid = topo.get_feature().get_reconstruction_plate_id()
            except Exception:
                pid = None
        if pid == int(plate_id):
            return topo
    return None


def compute_boundary_metrics(model, plate_id: int, age_ma: float, anchor_plate_id: int = 0) -> BoundaryMetrics:
    """
    Resolve the topologies for a plate and compute boundary composition and lengths.

    Parameters
    ----------
    model : ReconModel | iPlateModel
        Model providing rotation file and topology layer paths (via get_topologies()).
    plate_id : int
        Plate ID to analyse.
    age_ma : float
        Geological time (Ma).

    Returns
    -------
    BoundaryMetrics
    """
    _ensure_pygplates()

    topo_files = model.get_topologies()
    rot_model = model.rotation_model()
    # Respect anchor plate id to resolve relative topology correctly
    topo_model = pygplates.TopologicalModel(topo_files, rot_model, anchor_plate_id=anchor_plate_id)

    rtopo = _resolve_topology_for_plate(topo_model, age_ma, plate_id)
    if rtopo is None:
        return BoundaryMetrics(
            plate_id=int(plate_id),
            age_ma=age_ma,
            total_length_km=0.0,
            subduction_length_km=0.0,
            ridge_length_km=0.0,
            transform_length_km=0.0,
            other_length_km=0.0,
            subduction_fraction=0.0,
            centroid_lat=None,
            centroid_lon=None,
        )

    total_km = 0.0
    sz_km = 0.0
    ridge_km = 0.0
    tr_km = 0.0
    other_km = 0.0

    boundary_midpoints: List[pygplates.PointOnSphere] = []
    boundary_weights: List[float] = []

    for segment in rtopo.get_boundary_sub_segments():
        feature = segment.get_feature()
        resolved_feature = segment.get_resolved_feature()
        geometry = resolved_feature.get_geometry()
        if geometry is None:
            continue

        segment_km = _polyline_length_km(geometry)
        total_km += segment_km

        ftype = feature.get_feature_type()
        if ftype == pygplates.FeatureType.gpml_subduction_zone:
            sz_km += segment_km
        elif ftype == pygplates.FeatureType.gpml_mid_ocean_ridge:
            ridge_km += segment_km
        elif ftype == pygplates.FeatureType.gpml_transform:
            tr_km += segment_km
        else:
            other_km += segment_km

        points = getattr(geometry, 'get_points', lambda: None)()
        mids = _polyline_midpoints(points) if points else []
        for mp in mids:
            boundary_midpoints.append(mp)
            boundary_weights.append(max(segment_km / max(len(mids), 1), 1e-6))

    frac = (sz_km / total_km) if total_km > 0 else 0.0

    # Approximate centroid as weighted mean of boundary midpoints on sphere
    centroid_xyz = (0.0, 0.0, 0.0)
    for mp, w in zip(boundary_midpoints, boundary_weights):
        centroid_xyz = _add(centroid_xyz, _scale(_point_to_xyz(mp), w))
    centroid_xyz = _normalize(centroid_xyz)
    centroid_lat, centroid_lon = _xyz_to_lat_lon(centroid_xyz)

    return BoundaryMetrics(
        plate_id=int(plate_id),
        age_ma=age_ma,
        total_length_km=total_km,
        subduction_length_km=sz_km,
        ridge_length_km=ridge_km,
        transform_length_km=tr_km,
        other_length_km=other_km,
        subduction_fraction=frac,
        centroid_lat=centroid_lat,
        centroid_lon=centroid_lon,
    )


def compute_euler_pole(
    model,
    plate_id: int,
    age_ma: float,
    age_span: float,
    *,
    anchor_plate_id: int = 0,
    velocity_coefficient: float = 1.0,
    enforce_towards_subduction: bool = True,
    mode: str = 'segments',  # 'segments' (least-squares) or 'centroid' (pure centroid alignment)
    invert_sign: bool = False,
    debug: bool = False,
) -> Dict[str, Optional[float]]:
    """
    Compute an Euler pole for a plate that points its motion towards the densest
    distribution of subduction zones.

    Angular velocity (deg/Myr) is computed as:
        omega_deg_per_ma = subduction_fraction * velocity_coefficient

    Total rotation angle over an interval is:
        total_rotation_angle_deg = omega_deg_per_ma * age_span

    Returns a dict containing:
      - euler_lat, euler_lon: Euler pole position in degrees (None if indeterminate)
      - omega_deg_per_ma: angular speed (deg/Myr)
      - total_rotation_angle_deg: omega_deg_per_ma * age_span
      - subduction_fraction: ratio in [0,1]
      - metrics: BoundaryMetrics (as dict)
    """
    _ensure_pygplates()

    metrics = compute_boundary_metrics(model, plate_id, age_ma, anchor_plate_id=anchor_plate_id)
    if metrics.total_length_km <= 0 or metrics.centroid_lat is None:
        return {
            'euler_lat': None,
            'euler_lon': None,
            'omega_deg_per_ma': 0.0,
            'total_rotation_angle_deg': 0.0,
            'subduction_fraction': 0.0,
            'metrics': metrics.__dict__,
        }

    topo_files = model.get_topologies()
    rot_model = model.rotation_model()
    topo_model = pygplates.TopologicalModel(topo_files, rot_model, anchor_plate_id=anchor_plate_id)
    rtopo = _resolve_topology_for_plate(topo_model, age_ma, plate_id)
    if rtopo is None:
        return {
            'euler_lat': None,
            'euler_lon': None,
            'omega_deg_per_ma': 0.0,
            'total_rotation_angle_deg': 0.0,
            'subduction_fraction': 0.0,
            'metrics': metrics.__dict__,
        }

    centroid_point = pygplates.PointOnSphere(metrics.centroid_lat, metrics.centroid_lon)
    x = _point_to_xyz(centroid_point)

    # Accumulate tangent directions towards subduction midpoints, weighted by segment length
    sum_tangent = (0.0, 0.0, 0.0)
    per_seg_tangents: List[Tuple[float, float, float]] = []
    # Accumulator for global least-squares axis fit
    r_sum = (0.0, 0.0, 0.0)
    total_w = 0.0
    has_subduction = False

    for segment in rtopo.get_boundary_sub_segments():
        feature = segment.get_feature()
        if feature.get_feature_type() != pygplates.FeatureType.gpml_subduction_zone:
            continue
        resolved_feature = segment.get_resolved_feature()
        geometry = resolved_feature.get_geometry()
        if geometry is None:
            continue
        seg_len_km = _polyline_length_km(geometry)
        points = getattr(geometry, 'get_points', lambda: None)()
        mids = _polyline_midpoints(points) if points else []
        if not mids:
            continue
        has_subduction = True
        # Use mean of segment midpoints as direction target for this segment
        seg_mid_xyz = (0.0, 0.0, 0.0)
        for mp in mids:
            seg_mid_xyz = _add(seg_mid_xyz, _point_to_xyz(mp))
        seg_mid_xyz = _normalize(seg_mid_xyz)

        # Tangent direction at centroid pointing toward seg_mid (project seg_mid onto tangent plane at x)
        proj = _add(seg_mid_xyz, _scale(x, -_dot(seg_mid_xyz, x)))
        t = _normalize(proj)
        per_seg_tangents.append(t)
        w = max(seg_len_km, 1e-6)
        sum_tangent = _add(sum_tangent, _scale(t, w))

        # Build least-squares axis contribution at the subduction location:
        # Desired local direction at segment is projection of centroid-level desired direction onto tangent at seg.
        # This encourages velocities at subduction to follow the global desired motion.
        # We'll compute this once u_dir is available; temporarily store seg_mid and weight.
        # For now, accumulate placeholders; will add to r_sum after u_dir known.
        # Defer by storing in a list
        if 'seg_list' not in locals():
            seg_list = []
        seg_list.append((seg_mid_xyz, w))

    if not has_subduction:
        return {
            'euler_lat': None,
            'euler_lon': None,
            'omega_deg_per_ma': 0.0,
            'total_rotation_angle_deg': 0.0,
            'subduction_fraction': metrics.subduction_fraction,
            'metrics': metrics.__dict__,
        }

    t_avg = _normalize(sum_tangent)
    if t_avg == (0.0, 0.0, 0.0):
        return {
            'euler_lat': None,
            'euler_lon': None,
            'omega_deg_per_ma': 0.0,
            'total_rotation_angle_deg': 0.0,
            'subduction_fraction': metrics.subduction_fraction,
            'metrics': metrics.__dict__,
        }

    # Optionally ensure average motion points toward subduction on average
    if enforce_towards_subduction and per_seg_tangents:
        # Compute mean dot product with per-segment directions
        mean_dot = 0.0
        for t_seg in per_seg_tangents:
            mean_dot += (t_avg[0]*t_seg[0] + t_avg[1]*t_seg[1] + t_avg[2]*t_seg[2])
        mean_dot /= max(len(per_seg_tangents), 1)
        if mean_dot < 0.0:
            t_avg = (-t_avg[0], -t_avg[1], -t_avg[2])

    omega_dir = None
    if mode == 'centroid':
        # Force centroid velocity to aim at dominant subduction direction
        omega_dir = _normalize(_cross(x, t_avg))
    else:
        # Local least-squares fit using ONLY subduction segments (plate-local):
        # For each subduction midpoint r (unit), desired v at r points from centroid toward r along the surface.
        try:
            if 'seg_list' in locals() and seg_list:
                M = np.zeros((3, 3), dtype=float)
                b_vec = np.zeros(3, dtype=float)
                for r_xyz, w in seg_list:
                    # Cross-product matrix [r]_x: [r]_x v = r × v ; so A omega = omega × r = -[r]_x omega
                    rx, ry, rz = r_xyz
                    r_x = np.array([[0.0, -rz,  ry],
                                    [ rz, 0.0, -rx],
                                    [-ry,  rx, 0.0]], dtype=float)
                    # Desired direction at r: from centroid x to boundary r
                    d_raw = _add(x, _scale(r_xyz, -_dot(x, r_xyz)))
                    d = _normalize(d_raw)
                    d = (-d[0], -d[1], -d[2])
                    d_vec = np.array(d, dtype=float)

                    A = -r_x
                    M += w * (A.T @ A)
                    b_vec += w * (A.T @ d_vec)

                # Solve M omega = b
                omega = np.linalg.solve(M, b_vec)
                norm = np.linalg.norm(omega)
                if norm > 0:
                    omega_dir = tuple((omega / norm).tolist())
        except Exception:
            omega_dir = None

        if omega_dir is None:
            # Fallback to centroid-only axis
            omega_dir = _normalize(_cross(x, t_avg))

    # Determine sign so that instantaneous velocity v_dir = (omega_dir × x) aligns with t_avg
    v_dir = _normalize(_cross(omega_dir, x))
    align_sign = 1.0
    if v_dir != (0.0, 0.0, 0.0) and _dot(v_dir, t_avg) < 0.0:
        align_sign = -1.0

    # Choose sign by maximizing alignment at subduction midpoints with finite rotation over age_span
    # Evaluate both signs and pick the better one
    def alignment_score(sign: float) -> float:
        if 'seg_list' not in locals() or not seg_list:
            return 0.0
        score = 0.0
        angle_rad = (max(0.0, min(1.0, metrics.subduction_fraction)) * float(velocity_coefficient) * float(age_span) * sign) * math.pi / 180.0
        for r_xyz, w in seg_list:
            # Finite rotation of r about omega_dir by angle_rad
            r_rot = _rotate_about_axis(r_xyz, omega_dir, angle_rad)
            # Tangent direction of displacement at r ~ (r_rot - r) projected to tangent plane at r
            disp = (r_rot[0] - r_xyz[0], r_rot[1] - r_xyz[1], r_rot[2] - r_xyz[2])
            disp = _add(disp, _scale(r_xyz, -_dot(disp, r_xyz)))
            disp = _normalize(disp)
            # Desired direction at r: from centroid x toward r (outward)
            d_raw = _add(x, _scale(r_xyz, -_dot(x, r_xyz)))
            d = _normalize(d_raw)
            d = (-d[0], -d[1], -d[2])
            score += w * _dot(disp, d)
        return score

    score_pos = alignment_score(+1.0)
    score_neg = alignment_score(-1.0)
    align_sign = 1.0 if score_pos >= score_neg else -1.0

    e_lat, e_lon = _xyz_to_lat_lon(omega_dir)

    # Signed angular speed and total angle
    omega_deg_per_ma = align_sign * max(0.0, min(1.0, metrics.subduction_fraction)) * float(velocity_coefficient)
    if invert_sign:
        omega_deg_per_ma = -omega_deg_per_ma
    total_rotation_angle_deg = omega_deg_per_ma * float(age_span)

    if debug:
        def tup_fmt(v):
            return f"({v[0]:+.4f}, {v[1]:+.4f}, {v[2]:+.4f})"
        num_t = len(per_seg_tangents)
        v_dir = _normalize(_cross(omega_dir, x))
        mean_dot = 0.0
        if per_seg_tangents:
            for t_seg in per_seg_tangents:
                mean_dot += _dot(t_avg, t_seg)
            mean_dot /= max(len(per_seg_tangents), 1)
        print(f"[euler.debug] age={age_ma}Ma pid={plate_id} sub_frac={metrics.subduction_fraction:.3f} n_sub={num_t}")
        print(f"[euler.debug] centroid(lat,lon)=({metrics.centroid_lat:.3f},{metrics.centroid_lon:.3f}) x={tup_fmt(x)}")
        print(f"[euler.debug] t_avg={tup_fmt(t_avg)} mean_dot={mean_dot:+.4f} align_sign={align_sign:+.0f}")
        print(f"[euler.debug] pole(lat,lon)=({e_lat:.3f},{e_lon:.3f}) v_dir={tup_fmt(v_dir)} omega={omega_deg_per_ma:+.4f} deg/Myr invert_sign={invert_sign}")

    return {
        'euler_lat': e_lat,
        'euler_lon': e_lon,
        'omega_deg_per_ma': omega_deg_per_ma,
        'total_rotation_angle_deg': total_rotation_angle_deg,
        'subduction_fraction': metrics.subduction_fraction,
        'metrics': metrics.__dict__,
    }


def format_euler_rot_line(
    model,
    plate_id: int,
    age_ma: float,
    age_span: float,
    fixed_plate_id: int = 0,
    *,
    anchor_plate_id: int = 0,
    velocity_coefficient: float = 1.0,
    fmt_age: str = ".1f",
    fmt_val: str = ".4f",
    comment: str = "euler-auto",
    mode: str = 'segments',
    invert_sign: bool = False,
) -> Optional[str]:
    """
    Create a single .rot-format line using the computed Euler pole and rotation angle.

    Returns a string like:
        "<moving_pid> <age>  <lat>  <lon>  <angle>  <fixed_pid>  ! <comment>"

    Returns None if an Euler pole cannot be determined.
    """
    pole = compute_euler_pole(
        model,
        plate_id=plate_id,
        age_ma=age_ma,
        age_span=age_span,
        anchor_plate_id=anchor_plate_id,
        velocity_coefficient=velocity_coefficient,
        mode=mode,
        invert_sign=invert_sign,
    )

    e_lat = pole.get('euler_lat')
    e_lon = pole.get('euler_lon')
    angle = pole.get('total_rotation_angle_deg')

    if e_lat is None or e_lon is None or angle is None:
        return None

    age_txt = format(age_ma, fmt_age)
    lat_txt = format(float(e_lat), fmt_val)
    lon_txt = format(float(e_lon), fmt_val)
    ang_txt = format(float(angle), fmt_val)
    line = f"{int(plate_id)} {age_txt}  {lat_txt}  {lon_txt}  {ang_txt}  {int(fixed_plate_id):03d}  ! {comment}"
    return line


def print_euler_rot_line(
    model,
    plate_id: int,
    age_ma: float,
    age_span: float,
    fixed_plate_id: int = 0,
    *,
    anchor_plate_id: int = 0,
    velocity_coefficient: float = 1.0,
    fmt_age: str = ".1f",
    fmt_val: str = ".4f",
    comment: str = "euler-auto",
    debug: bool = False,
    mode: str = 'segments',
    invert_sign: bool = False,
) -> Optional[str]:
    """
    Print and return a .rot-format line for copy/paste.
    Returns None if an Euler pole cannot be determined.
    """
    # Optional debug invoke of compute_euler_pole
    if debug:
        _ = compute_euler_pole(
            model,
            plate_id,
            age_ma,
            age_span,
            anchor_plate_id=anchor_plate_id,
            velocity_coefficient=velocity_coefficient,
            enforce_towards_subduction=True,
            debug=True,
            mode=mode,
            invert_sign=invert_sign,
        )

    line = format_euler_rot_line(
        model,
        plate_id,
        age_ma,
        age_span,
        fixed_plate_id,
        anchor_plate_id=anchor_plate_id,
        velocity_coefficient=velocity_coefficient,
        fmt_age=fmt_age,
        fmt_val=fmt_val,
        comment=comment,
        mode=mode,
        invert_sign=invert_sign,
    )
    if line is not None:
        print(line)
    # return line


