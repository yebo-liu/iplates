import os
try:
    import pygplates
except ImportError:
    pygplates = None

try:
    import geopandas as gpd
except ImportError:
    gpd = None


def resolve_topos(topo_model,reconstruction_time,write_path):
    '''
    A very crude function to generate the boundaries.
    '''
    # Ensure output directory exists
    if write_path and not os.path.isdir(write_path):
        os.makedirs(write_path, exist_ok=True)

    l_sz_features = []
    r_sz_features = []
    np_sz_features = []
    ridge_features = []
    tr_features = []
    other_features = []
    topo_sections = topo_model.topological_snapshot(reconstruction_time).get_resolved_topologies()
    for topo in topo_sections:
        for segment in topo.get_boundary_sub_segments():
            if segment.get_feature().get_feature_type() == pygplates.FeatureType.gpml_mid_ocean_ridge:
                    ridge_features.append(segment.get_resolved_feature())
            elif segment.get_feature().get_feature_type() == pygplates.FeatureType.gpml_transform:
                tr_features.append(segment.get_resolved_feature())
            elif segment.get_feature().get_feature_type() == pygplates.FeatureType.gpml_subduction_zone:
                polarity_property = segment.get_feature().get(
                    pygplates.PropertyName.create_gpml('subductionPolarity'))
                if polarity_property:
                    polarity = polarity_property.get_value().get_content()
                    if polarity == 'Left':
                        l_sz_features.append(segment.get_resolved_feature())
                    elif polarity == 'Right':
                        r_sz_features.append(segment.get_resolved_feature())
                    elif polarity == 'Unknown':
                        np_sz_features.append(segment.get_resolved_feature())
                else:
                    # No polarity property; treat as unknown
                    np_sz_features.append(segment.get_resolved_feature())
                    print('Encountered subduction segment without polarity; classified as unknown.')
            else:
                other_features.append(segment.get_resolved_feature())
    pygplates.FeatureCollection(l_sz_features).write(f'{write_path}/sz_L_{reconstruction_time}.shp')
    pygplates.FeatureCollection(r_sz_features).write(f'{write_path}/sz_R_{reconstruction_time}.shp')
    pygplates.FeatureCollection(np_sz_features).write(f'{write_path}/sz_np_{reconstruction_time}.shp')
    pygplates.FeatureCollection(ridge_features).write(f'{write_path}/ridge_{reconstruction_time}.shp')
    pygplates.FeatureCollection(tr_features).write(f'{write_path}/tr_{reconstruction_time}.shp')
    pygplates.FeatureCollection(other_features).write(f'{write_path}/other_{reconstruction_time}.shp')

def plot_boundaries(reconstruction_time,input_path):
    sz_R_path = f'{input_path}/sz_R_{reconstruction_time}.shp'
    if os.path.isfile(sz_R_path):
        gdf_r = gpd.read_file(sz_R_path)
        fig.plot(gdf_r, style='f30p/4p+r+t+o7p', fill='#f26c64', pen='1p,#f26c64')
    
    sz_L_path = f'{input_path}/sz_L_{reconstruction_time}.shp'
    if os.path.isfile(sz_L_path):
        gdf_l = gpd.read_file(sz_L_path)
        fig.plot(gdf_l, style='f30p/4p+l+t+o7p', fill='#f26c64', pen='1p,#f26c64')

    sz_np_path = f'{input_path}/sz_np_{reconstruction_time}.shp'
    if os.path.isfile(sz_np_path):
        gdf_np = gpd.read_file(sz_np_path)
        fig.plot(gdf_np,  pen='1p,#f26c64')
        
        
    mor_path = f'{input_path}/ridge_{reconstruction_time}.shp'
    if os.path.isfile(mor_path):
        gdf_mor = gpd.read_file(mor_path)
        fig.plot(data=gdf_mor, pen='1p,#028100')
        
    tr_path = f'{input_path}/tr_{reconstruction_time}.shp'    
    if os.path.isfile(tr_path):
        gdf_tr = gpd.read_file(tr_path)
        fig.plot(data=gdf_tr, pen='1p,#d1a98c')

    other_path = f'{input_path}/other_{reconstruction_time}.shp'
    if os.path.isfile(other_path):
        gdf_other = gpd.read_file(other_path)
        fig.plot(data=gdf_other, pen='1p,#028100')
