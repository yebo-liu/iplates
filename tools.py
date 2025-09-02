import os
import pygplates
import geopandas as gpd
import pandas as pd
from pathlib import Path
from typing import Union, List
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon

def feature_collection_to_gdf(feature_collection, plate_id_col='PLATEID1'):
    """
    Convert a pygplates FeatureCollection to a GeoPandas GeoDataFrame.
    
    Parameters:
    -----------
    feature_collection : pygplates.FeatureCollection
        The feature collection to convert
    plate_id_col : str, optional
        Name for the plate ID column (default: 'PLATEID1')
        
    Returns:
    --------
    geopandas.GeoDataFrame
        A GeoDataFrame with all features and their properties
    """
    features_data = []
    
    for feature in feature_collection:
        # Get geometry
        geometry = feature.get_geometry()
        if geometry is None:
            continue
            
        # Convert pygplates geometry to shapely
        shapely_geom = pygplates_to_shapely_geometry(geometry)
        if shapely_geom is None:
            continue
            
        # Get reconstruction plate ID
        plate_id = feature.get_reconstruction_plate_id()
        if plate_id is None:
            plate_id = 0
            
        # Get feature properties
        properties = {}
        properties[plate_id_col] = plate_id
        
        # Get shapefile attributes
        for attr_name in feature.get_shapefile_attribute_names():
            try:
                attr_value = feature.get_shapefile_attribute(attr_name)
                if attr_value is not None:
                    properties[attr_name] = attr_value
            except:
                continue
                
        # Get GPML properties
        for prop_name in feature.get_property_names():
            try:
                prop = feature.get(prop_name)
                if prop:
                    prop_value = prop.get_value()
                    if hasattr(prop_value, 'get_content'):
                        properties[str(prop_name)] = prop_value.get_content()
                    else:
                        properties[str(prop_name)] = str(prop_value)
            except:
                continue
                
        # Get feature name
        feature_name = feature.get_name()
        if feature_name:
            properties['NAME'] = feature_name
            
        # Get feature type
        feature_type = feature.get_feature_type()
        if feature_type:
            properties['FEATURE_TYPE'] = str(feature_type)
            
        # Add geometry and properties to our data
        feature_data = properties.copy()
        feature_data['geometry'] = shapely_geom
        features_data.append(feature_data)
    
    if not features_data:
        # Return empty GeoDataFrame with proper structure
        return gpd.GeoDataFrame(columns=[plate_id_col, 'geometry'])
        
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(features_data)
    
    # Set CRS to WGS84 (standard for plate tectonic data)
    gdf.set_crs(epsg=4326, inplace=True)
    
    return gdf

def pygplates_to_shapely_geometry(geometry):
    """
    Convert a pygplates geometry to a shapely geometry.
    """
    if isinstance(geometry, pygplates.PointOnSphere):
        lat, lon = geometry.to_lat_lon()
        return Point(lon, lat)  # shapely uses (x, y) = (lon, lat)
        
    elif isinstance(geometry, pygplates.MultiPointOnSphere):
        points = []
        for pt in geometry:
            lat, lon = pt.to_lat_lon()
            points.append(Point(lon, lat))
        return MultiPoint(points)
        
    elif isinstance(geometry, pygplates.PolylineOnSphere):
        coords = []
        for pt in geometry:
            lat, lon = pt.to_lat_lon()
            coords.append((lon, lat))
        return LineString(coords)
        
    elif isinstance(geometry, pygplates.MultiPolylineOnSphere):
        lines = []
        for line in geometry:
            coords = []
            for pt in line:
                lat, lon = pt.to_lat_lon()
                coords.append((lon, lat))
            lines.append(LineString(coords))
        return MultiLineString(lines)
        
    elif isinstance(geometry, pygplates.PolygonOnSphere):
        coords = []
        for pt in geometry:
            lat, lon = pt.to_lat_lon()
            coords.append((lon, lat))
        return Polygon(coords)
        
    elif isinstance(geometry, pygplates.MultiPolygonOnSphere):
        polygons = []
        for poly in geometry:
            coords = []
            for pt in poly:
                lat, lon = pt.to_lat_lon()
                coords.append((lon, lat))
            polygons.append(Polygon(coords))
        return MultiPolygon(polygons)
        
    else:
        return None

def shapely_to_pygplates_geometry(geometry):
    """
    Convert a shapely geometry to a pygplates geometry.
    """
    if isinstance(geometry, Point):
        return pygplates.PointOnSphere(geometry.y, geometry.x)  # lat, lon
    elif isinstance(geometry, MultiPoint):
        return pygplates.MultiPointOnSphere([pygplates.PointOnSphere(pt.y, pt.x) for pt in geometry.geoms])
    elif isinstance(geometry, LineString):
        return pygplates.PolylineOnSphere([(pt[1], pt[0]) for pt in geometry.coords])  # (lat, lon)
    elif isinstance(geometry, MultiLineString):
        return pygplates.MultiPolylineOnSphere([
            pygplates.PolylineOnSphere([(pt[1], pt[0]) for pt in line.coords])
            for line in geometry.geoms
        ])
    elif isinstance(geometry, Polygon):
        # Only use exterior ring for now
        return pygplates.PolygonOnSphere([(pt[1], pt[0]) for pt in geometry.exterior.coords])
    elif isinstance(geometry, MultiPolygon):
        return pygplates.MultiPolygonOnSphere([
            pygplates.PolygonOnSphere([(pt[1], pt[0]) for pt in poly.exterior.coords])
            for poly in geometry.geoms
        ])
    else:
        return None

def gdf_to_feature_collection(
    gdf, 
    feature_type='gpml:UnclassifiedFeature', 
    reconstruction_plate_id_col='PLATEID1'
):
    features = []
    for idx, row in gdf.iterrows():
        geometry = shapely_to_pygplates_geometry(row.geometry)
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
    # print(f"File: {rot_file_path}")
    print(f"Active interval: {oldest['age_ma']} Ma -> {youngest['age_ma']} Ma   ")
    print("")

    print("-- Rotation Lines --")
    for occ in occurrences:
        print(f"{occ['raw']} ==== [Line {occ['line_num']}]")

    print("")
    chain_info = get_plate_chain(rot_file_path, plate_id, age)
    if chain_info is None:
        age_text = f" at {age} Ma" if age is not None else ""
        print(f"-- Fixed-plate chain{age_text} --")
        print("  Not available.")
        return

    chain_repr = " -> ".join(f"{pid:03d}" for pid in chain_info['chain'])
    print(f"-- Fixed-plate chain at {chain_info['age']} Ma --")
    print(f"  {chain_repr}")

