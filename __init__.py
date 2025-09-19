__version__ = "0.1.0"
VERSION = __version__
AUTHOR = "Yebo Liu"

# Main classes
from .loader import ReconModel, iPlateModel, loadW
# from .reconstruction import *

# Import modules to make them accessible
from . import loader
from . import reconstruction
from . import tools
from . import plot
from . import euler

# Key functions from tools
from .tools import (
    print_plate_chain,
    get_plate_active_interval,
    get_plate_chain,
    find_plate_occurrences_in_rot_file,
    get_pids_fixed_to_target,
    feature_collection_to_gdf,
    gdf_to_feature_collection
)

# Euler exports
from .euler import (
    BoundaryMetrics,
    compute_boundary_metrics,
    compute_euler_pole,
    format_euler_rot_line,
    print_euler_rot_line,
)

# Make key functions available at package level
__all__ = [
    'ReconModel',
    'iPlateModel', 
    'loadW',
    'print_plate_chain',
    'get_plate_active_interval',
    'get_plate_chain',
    'find_plate_occurrences_in_rot_file',
    'get_pids_fixed_to_target',
    'feature_collection_to_gdf',
    'gdf_to_feature_collection'
    ,
    'BoundaryMetrics',
    'compute_boundary_metrics',
    'compute_euler_pole',
    'format_euler_rot_line',
    'print_euler_rot_line'
]