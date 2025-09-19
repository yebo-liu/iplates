import os
import pygplates
import gplately
import geopandas as gpd
import pandas as pd
from pathlib import Path
from typing import Union, List
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon

def _normalize_feature_collection_input(feature_collection):
    """
    Normalize different input types to a list of pygplates.Feature objects.
    
    Handles:
    1. str: path to a file -> load as FeatureCollection
    2. list of str: list of paths -> load as FeatureCollection (pygplates supports this)
    3. list of tuples: [(Feature, reconstructed_geoms), ...] -> extract Features and set geometries
    4. Already a FeatureCollection or list of Features -> return as-is
    """
    
    # Type 1: String path
    if isinstance(feature_collection, str):
        try:
            return pygplates.FeatureCollection(feature_collection)
        except Exception as e:
            raise ValueError(f"Cannot load file '{feature_collection}': {e}")
    
    # Type 2: List of strings (file paths) - pygplates supports this directly
    elif isinstance(feature_collection, list) and len(feature_collection) > 0:
        if all(isinstance(item, str) for item in feature_collection):
            # All strings - treat as file paths, use first one
            try:
                return pygplates.FeatureCollection(feature_collection[0])
            except Exception as e:
                raise ValueError(f"Cannot load file '{feature_collection[0]}': {e}")
        
        # Type 3: List of tuples [(Feature, reconstructed_geoms), ...]
        # The code below extract the reconstructed geometries and set them to the features.
        elif all(isinstance(item, (tuple, list)) and len(item) >= 2 for item in feature_collection):
            features = []
            for feature, recon_geoms in feature_collection:
                # Set reconstructed geometry on the feature
                if recon_geoms and hasattr(recon_geoms[0], 'get_reconstructed_geometry'):
                    try:
                        real_reconstructed_geometry = [geom.get_reconstructed_geometry() for geom in recon_geoms]
                        feature.set_geometry(real_reconstructed_geometry)
                    except Exception:
                        pass  # Skip if geometry setting fails
                features.append(feature)
            return features
        
        # Type 4: Already a list of Features
        elif all(hasattr(item, 'get_geometry') for item in feature_collection):
            return feature_collection
    
    # Type 5: Single FeatureCollection object
    elif hasattr(feature_collection, '__iter__') and hasattr(feature_collection, '__getitem__'):
        try:
            return feature_collection
        except Exception:
            pass
    
    # Type 6: Single Feature
    elif hasattr(feature_collection, 'get_geometry'):
        return [feature_collection]
    
    raise ValueError(f"Unsupported input type: {type(feature_collection)}. "
                   f"Expected: str, list of str, list of tuples, FeatureCollection, or list of Features")

def feature_collection_to_gdf(feature_collection):
    """
    Convert a pygplates FeatureCollection to a GeoPandas GeoDataFrame.
    
    Parameters:
    -----------
    feature_collection can be:
        1. a path to a file that can be read by pygplates
        2. the output of a pygplates.ReconstructSnapshot.get_reconstructed_features()

    based on the type of feature_collection, the code will handle the data differently.
    Returns:
    --------
    geopandas.GeoDataFrame
        A GeoDataFrame with all features and their properties
    """
    feature_collection = _normalize_feature_collection_input(feature_collection)

    features_data = []
    for feature in feature_collection:

        # Get and set feature properties
        properties = {}
        
        # Get shapefile attributes (some gpmls doesn't have these attributes)
        shapefile_attributes = feature.get_shapefile_attributes()
        if shapefile_attributes:
            # shapefile_attributes is a dict, so iterate through it
            for attr_name, attr_value in shapefile_attributes.items():
                properties[attr_name] = attr_value

        # when there are no shapefile attributes, we will try to get the gpml properties
        # TODO: for now, the code only extracts a few essential gpml properties, should extract all and 
        # map their names to the conventional ones.
        else:
            # continue
            reconstructed_pid = feature.get_reconstruction_plate_id()
            feature_id = feature.get_feature_id()
            feature_type = feature.get_feature_type()
            from_age = feature.get_valid_time()[0]
            to_age = feature.get_valid_time()[1]
            feature_name = feature.get_name()
            geom_import_age = feature.get_geometry_import_time()

            properties['PLATEID1'] = reconstructed_pid
            properties['FEATURE_ID'] = feature_id
            properties['GPGIM_TYPE'] = feature_type
            properties['FROMAGE'] = from_age
            properties['TOAGE'] = to_age
            properties['NAME'] = feature_name
            properties['IMPORT_AGE'] = geom_import_age

        # Get geometry
        geometry = feature.get_geometry()
        if geometry is None:
            continue
            
        # Convert pygplates geometry to shapely
        shapely_geom = gplately.geometry.pygplates_to_shapely(geometry)
        if shapely_geom is None:
            continue

        properties['geometry'] = shapely_geom
        # Add geometry and properties to our data
        features_data.append(properties)
        
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(features_data)
    
    # Set CRS to WGS84 (standard for plate tectonic data)
    gdf.set_crs(epsg=4326, inplace=True)
    
    return gdf


def gdf_to_feature_collection(
    gdf, 
    feature_type='gpml:UnclassifiedFeature', 
    reconstruction_plate_id_col='PLATEID1'
):
    features = []
    for idx, row in gdf.iterrows():
        geometry = gplately.geometry.shapely_to_pygplates(row.geometry)
        if geometry is None:
            print(f"Warning: Unsupported geometry at index {idx}, skipping.")
            continue

        # Correct feature creation
        if feature_type == 'gpml:UnclassifiedFeature':
            feature = pygplates.Feature()
        else:
            feature_type_obj = pygplates.FeatureType.create_gpml(feature_type.split(':')[-1])
            feature = pygplates.Feature(feature_type_obj)

        feature.set_geometry(geometry)

        # Set reconstruction plate ID if available
        if reconstruction_plate_id_col in row:
            try:
                plate_id = int(row[reconstruction_plate_id_col])
            except Exception:
                plate_id = 0
            feature.set_reconstruction_plate_id(plate_id)

        # Add properties as feature properties
        # TODO: this doesn't really work, perhaps because set_shapefile_attribute doesn't work with gpml?
        for col in gdf.columns:
            if col not in ['geometry', reconstruction_plate_id_col]:
                value = row[col]
                if pd.notnull(value):
                    feature.set_shapefile_attribute(str(col), str(value))

        features.append(feature)

    return pygplates.FeatureCollection(features)

def traverse_sub_tree(edge, list_to_append=None):

    if list_to_append is None:
        list_to_append = []

    for child_edge in edge.get_child_edges():
        list_to_append.append(child_edge.get_moving_plate_id())
        traverse_sub_tree(child_edge, list_to_append)

    return list_to_append


def get_pids_fixed_to_target(rot_file, age, target_pid=None):
    '''
    This function returns the list of plate IDs that are fixed to the target plate at a given age.
    '''
    rot_model = pygplates.RotationModel(rot_file)
    for edge in rot_model.get_reconstruction_tree(age).get_edges():
        if edge.get_moving_plate_id() == target_pid:
            return traverse_sub_tree(edge)


def _parse_rot_line(line_text: str):
    """
    Parse a single .rot line into its components.

    Expected layout (whitespace separated, optional trailing comment prefixed by '!'):
        <moving_plate_id> <age_Ma> <euler_lat> <euler_lon> <angle_deg> <fixed_plate_id> !<comment>

    Returns a dict with parsed fields or None if the line does not conform.
    """
    # Remove newline and keep raw copy
    raw = line_text.rstrip("\n")

    # Separate optional comment (anything after '!')
    excl_index = raw.find('!')
    if excl_index >= 0:
        head = raw[:excl_index].strip()
        comment = raw[excl_index + 1 :].strip()
    else:
        head = raw.strip()
        comment = ''

    if not head:
        return None

    parts = head.split()
    # We need at least 6 columns
    if len(parts) < 6:
        return None

    try:
        moving_plate_id = int(parts[0])
        age_ma = float(parts[1])
        euler_lat = float(parts[2])
        euler_lon = float(parts[3])
        angle_deg = float(parts[4])
        fixed_plate_id = int(parts[5])
    except Exception:
        return None

    return {
        'moving_plate_id': moving_plate_id,
        'age_ma': age_ma,
        'euler_lat': euler_lat,
        'euler_lon': euler_lon,
        'angle_deg': angle_deg,
        'fixed_plate_id': fixed_plate_id,
        'comment': comment,
        'raw': raw,
    }


def find_plate_occurrences_in_rot_file(rot_file_path: str, plate_id: int):
    """
    Scan a .rot file for all lines belonging to a moving PLATEID.

    Parameters
    ----------
    rot_file_path : str
        Absolute or relative path to the .rot file
    plate_id : int
        Moving plate ID of interest

    Returns
    -------
    list[dict]
        List of occurrences sorted by age (ascending), each dict contains:
        { 'line_num', 'raw', 'age_ma', 'fixed_plate_id', ...all parsed fields }
    """
    occurrences = []
    with open(rot_file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line_index, line_text in enumerate(f, start=1):
            parsed = _parse_rot_line(line_text)
            if not parsed:
                continue
            if parsed['moving_plate_id'] == int(plate_id):
                parsed_with_line = parsed.copy()
                parsed_with_line['line_num'] = line_index
                occurrences.append(parsed_with_line)

    # Sort youngest->oldest (age ascending)
    occurrences.sort(key=lambda item: item['age_ma'])
    return occurrences


def get_plate_active_interval(rot_file_path: str, plate_id: int):
    """
    Determine the youngest and oldest ages where the given PLATEID appears in a .rot file.

    Returns
    -------
    dict | None
        None if PLATEID not found. Otherwise a dict with:
        {
            'youngest_age': float,
            'oldest_age': float,
            'youngest_occurrence': {'line_num': int, 'raw': str, ...},
            'oldest_occurrence': {'line_num': int, 'raw': str, ...},
            'all_occurrences': list[dict]
        }
    """
    occurrences = find_plate_occurrences_in_rot_file(rot_file_path, plate_id)
    if not occurrences:
        return None

    youngest = occurrences[0]
    oldest = occurrences[-1]
    return {
        'youngest_age': youngest['age_ma'],
        'oldest_age': oldest['age_ma'],
        'youngest_occurrence': youngest,
        'oldest_occurrence': oldest,
        'all_occurrences': occurrences,
    }


def _infer_chain_age_for_plate(rot_file_path: str, plate_id: int) -> float | None:
    """
    Prefer an age close to present day if available, otherwise choose the youngest
    age at which the plate appears in the .rot file.
    """
    occurrences = find_plate_occurrences_in_rot_file(rot_file_path, plate_id)
    if not occurrences:
        return None
    # Try exact 0.0 if present
    for occ in occurrences:
        if abs(occ['age_ma'] - 0.0) < 1e-9:
            return 0.0
    # Fallback to youngest available
    return occurrences[0]['age_ma']


def get_plate_chain(rot_file_path: str, plate_id: int, age: float | None = None):
    """
    Build the fixed-plate chain for a PLATEID at a given age using pygplates' reconstruction tree.

    Parameters
    ----------
    rot_file_path : str
        Path to .rot file
    plate_id : int
        Moving plate ID of interest
    age : float | None
        Geological age (Ma). If None, use 0.0 if available, else the plate's youngest age.

    Returns
    -------
    dict | None
        None if PLATEID not present at the chosen age or tree cannot be built.
        Otherwise a dict containing:
        {
          'age': float,
          'chain': list[int],            # eg [60302, 603, 616, 801, 000]
          'edges': list[pygplates.ReconstructionTreeEdge],
        }
    """
    chosen_age = age if age is not None else _infer_chain_age_for_plate(rot_file_path, plate_id)
    if chosen_age is None:
        return None

    rot_model = pygplates.RotationModel(rot_file_path)
    tree = rot_model.get_reconstruction_tree(chosen_age)

    # Find the edge corresponding to our plate_id
    start_edge = None
    for edge in tree.get_edges():
        if edge.get_moving_plate_id() == int(plate_id):
            start_edge = edge
            break

    if start_edge is None:
        return None

    chain_ids = [start_edge.get_moving_plate_id()]
    chain_edges = [start_edge]

    # Walk up the tree to the anchor (root)
    current_edge = start_edge
    while True:
        parent_edge = current_edge.get_parent_edge()
        if parent_edge is None:
            # Include the fixed plate id at this edge to close the chain
            # The fixed plate of the topmost edge is typically the anchor plate
            anchor_id = current_edge.get_fixed_plate_id()
            chain_ids.append(anchor_id)
            break
        chain_ids.append(parent_edge.get_moving_plate_id())
        chain_edges.append(parent_edge)
        current_edge = parent_edge

    return {'age': chosen_age, 'chain': chain_ids, 'edges': chain_edges}


def print_plate_chain(rot_file_path: str, plate_id: int, age: float | None = None):
    """
    Print a detailed report for a PLATEID including:
    - Active age interval
    - All .rot lines for that PLATEID (with line numbers)
    - Fixed-plate chain at a chosen (or inferred) age
    """
    occurrences = find_plate_occurrences_in_rot_file(rot_file_path, plate_id)
    if not occurrences:
        print(f"PLATEID {plate_id} not found in {rot_file_path}")
        return

    youngest = occurrences[0]
    oldest = occurrences[-1]

    print(f"=== Plate {plate_id} ===")
    print(f"Active interval: {oldest['age_ma']} Ma -> {youngest['age_ma']} Ma   ")
    print("")

    chain_info = get_plate_chain(rot_file_path, plate_id, age)
    if chain_info is None:
        age_text = f" at {age} Ma" if age is not None else ""
        print(f"-- Plate chain{age_text} --")
        print("  Not available.")
        return

    chain_repr = " -> ".join(f"{pid:03d}" for pid in chain_info['chain'])
    print(f"-- Plate chain at {chain_info['age']} Ma --")
    print(f"  {chain_repr}")
    print("")

    print("-- Rotation Lines --")
    print(f'Line Number: {occurrences[0]["line_num"]} - {occurrences[-1]["line_num"]}')
    for occ in occurrences:
        print(f"{occ['raw']}")

def reanchor_plate(
    rot_file_path: str,
    moving_plate_id: int,
    new_fixed_plate_id: int,
    ages: List[float] | None = None,
    wrap_angle_on_path_error: bool = False,
    fmt_age: str = ".1f",
    fmt_val: str = ".4f",
):
    """
    Recompute and PRINT rotation lines for a moving plate relative to a different fixed plate (anchor),
    following your previous workflow using `pygplates.RotationModel.get_rotation`.

    For each selected age, computes the finite rotation of `moving_plate_id` relative to
    `new_fixed_plate_id` and prints standardized lines. No files are written.

    Parameters
    ----------
    rot_file_path : str
        Path to source .rot file (only used to initialize the rotation model).
    moving_plate_id : int
        Plate whose reference frame should change.
    new_fixed_plate_id : int
        The new anchor (fixed plate) to reference against.
    ages : list[float] | None
        If provided, compute lines for these ages. If None, use all ages present for the plate
        in the file (via find_plate_occurrences_in_rot_file).
    wrap_angle_on_path_error : bool
        Apply the angle wrapping you used (±360 shift based on sign).
    fmt_age : str
        Format spec for the age field (default .1f).
    fmt_val : str
        Format spec for lat/lon/angle (default .4f).

    Returns
    -------
    list[str]
        The printed lines (also returned for convenience).
    """
    # Gather target ages
    if ages is None:
        occ = find_plate_occurrences_in_rot_file(rot_file_path, moving_plate_id)
        if not occ:
            return {
                'plate_id': int(moving_plate_id),
                'new_fixed_plate_id': int(new_fixed_plate_id),
                'num_target_ages': 0,
                'changed_lines': 0,
                'output_file': None,
                'preview': [],
            }
        ages = [o['age_ma'] for o in occ]

    # Load rotation model once
    rot_model = pygplates.RotationModel(rot_file_path)

    # Compute new lines and print
    lines_out: List[str] = []
    for age in ages:
        finite_rotation = rot_model.get_rotation(age, moving_plate_id, anchor_plate_id=new_fixed_plate_id)
        pole_lat, pole_lon, pole_angle = finite_rotation.get_lat_lon_euler_pole_and_angle_degrees()

        if wrap_angle_on_path_error:
            # Reflect the angle by ±360 like your previous snippet
            if pole_angle >= 0:
                pole_angle = pole_angle - 360.0
            else:
                pole_angle = pole_angle + 360.0

        # Format similar to your prints, keeping extra spaces between numeric columns
        # Example: "305 380.0  12.3456  123.4567  -45.6789  302  !"
        age_txt = format(age, fmt_age)
        lat_txt = format(pole_lat, fmt_val)
        lon_txt = format(pole_lon, fmt_val)
        ang_txt = format(pole_angle, fmt_val)
        line = f"{int(moving_plate_id)} {age_txt}  {lat_txt}  {lon_txt}  {ang_txt}  {int(new_fixed_plate_id):03d}  !"
        print(line)
        lines_out.append(line)

    return lines_out