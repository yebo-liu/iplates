import pygplates
from .tools import feature_collection_to_gdf
import numpy as np
import pandas as pd
import pygmt



def highlight_pids_on_map(model, pids,highlight_type='COBs',age = 0.0,anchor_plate_id = 0):
    '''This function do the following things:
    1. Reconstruct the given type of features (default COBs) at the given age
    2. Highlight the given pids on the map
    '''
    if highlight_type == 'COBs':
        reconstructable_features = model.get_COBs()
    elif highlight_type == 'Coastlines':
        reconstructable_features = model.get_coastlines()

    rot_file = model.get_rotation_model()
    model_recon = pygplates.ReconstructModel(reconstructable_features,rot_file,anchor_plate_id)
    recon_snapshot = model_recon.reconstruct_snapshot(age)
    reconstructed_features = recon_snapshot.get_reconstructed_features()
    gdf_recon = feature_collection_to_gdf(reconstructed_features)
    _plot_map_highlighting_pids(gdf_recon, pids)

def _plot_map_highlighting_pids(gdf_recon, pids):
    '''
    Plots gdf_recon. For geometries with given pids, the color will be brighter. The rest
    of geometries will be dull.
    '''
    if gdf_recon is None or len(gdf_recon) == 0:
        raise ValueError("gdf_recon is empty; nothing to plot")

    # Normalize pids to a set of ints for fast lookup
    try:
        pid_set = {int(pid) for pid in pids}
    except Exception:
        pid_set = set(pids)

    # Determine highlight mask; default to no highlight if column missing
    if 'PLATEID1' in gdf_recon.columns:
        try:
            mask = gdf_recon['PLATEID1'].astype(int).isin(pid_set)
        except Exception:
            mask = gdf_recon['PLATEID1'].isin(pid_set)
    else:
        mask = gdf_recon.assign(_tmp=False)['_tmp']

    base_fill = '#90A4AE'
    base_edge = '0.5p,#90A4AE'
    hi_edge = '1p,#D84315'
    hi_fill = '#FF7043'



    base_df = gdf_recon[~mask].copy()
    hi_df = gdf_recon[mask].copy()

    fig = pygmt.Figure()
 
    projection = 'W0/12c'  # Mollweide, central meridian 0, width 18 cm
    region = 'g'           # global
    pygmt.config(MAP_GRID_PEN_PRIMARY="0.6p,gray60,-")
    pygmt.config(MAP_GRID_PEN_SECONDARY="0.6p,gray60,-")
    pygmt.config(MAP_FRAME_PEN="0.6p,gray80")

    fig.basemap(region=region, projection=projection, frame='30g30')  

    if base_df is not None:
        fig.plot(data=base_df, pen=base_edge, fill=base_fill,transparency=50)

    if hi_df is not None:
        fig.plot(data=hi_df, pen=hi_edge, fill=hi_fill)

    fig.show(width=1000)