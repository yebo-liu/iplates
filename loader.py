import json
import os
from pathlib import Path
from typing import Dict, Optional, List, Union

import pygplates


class ReconModel:
    '''
    Minimal JSON-driven plate model loader (no downloading, local files only).

    JSON schema (minimal):
    {
      "BigTime": 18000,                   # optional, defaults None
      "SmallTime": 0,                      # optional, defaults None
      "Rotations": "/abs/path/model.rot",          # required (local path)
      "Layers": {                                     # optional named feature layers
        "StaticPolygons": "/abs/path/static.gpml",
        "Coastlines": "/abs/path/coastlines.gpml"
      },
      "TimeDepRasters": { }                # optional; stored as-is (paths)
    }
    '''

    def __init__(
        self,
        rotation_file: str,
        layers: Optional[Dict[str, object]] = None,
        big_time: Optional[float] = None,
        small_time: Optional[float] = None,
        time_dep_rasters: Optional[Dict[str, object]] = None,
        base_dir: Optional[str] = None,
    ) -> None:
        # Paths and metadata
        self._base_dir = str(Path(base_dir).resolve()) if base_dir else None
        self._rotation_file = _resolve_local(rotation_file, self._base_dir)
        self._layer_files: Dict[str, object] = layers or {}
        self._big_time = big_time
        self._small_time = small_time
        self._time_dep_rasters = time_dep_rasters or {}

        # Validation
        if not os.path.isfile(self._rotation_file):
            raise FileNotFoundError(f"Rotation file not found: {self._rotation_file}")

        # Validate layer file paths (strings or list of strings)
        for name, spec in self._layer_files.items():
            if isinstance(spec, str):
                _assert_local_file(spec, f"Layer '{name}'", self._base_dir)
            elif isinstance(spec, list):
                if not spec:
                    raise ValueError(f"Layer '{name}' is an empty list.")
                for p in spec:
                    if not isinstance(p, str):
                        raise ValueError(f"Layer '{name}' must be a string or list of strings.")
                    _assert_local_file(p, f"Layer '{name}' entry", self._base_dir)
            else:
                raise ValueError(f"Layer '{name}' must be a string path or list of string paths.")

        # Cache for rotation model only (layers are returned as paths)
        self._rotation_model: Optional[pygplates.RotationModel] = None

    @classmethod
    def from_json(cls, config_path: str) -> "ReconModel":
        cfg_path = str(Path(config_path).expanduser())
        with open(cfg_path, 'r', encoding='utf-8') as f:
            cfg = json.load(f)

        # Accept both camel-case and lower-case variants for robustness
        rotation_file = cfg.get('Rotations') or cfg.get('rotations') or cfg.get('rotation_file')
        if not rotation_file:
            raise ValueError("'Rotations' (path to .rot) is required in the config JSON.")

        return cls(
            rotation_file=rotation_file,
            layers=cfg.get('Layers') or cfg.get('layers') or {},
            big_time=cfg.get('BigTime') or cfg.get('big_time'),
            small_time=cfg.get('SmallTime') or cfg.get('small_time'),
            time_dep_rasters=cfg.get('TimeDepRasters') or cfg.get('time_dep_rasters') or {},
            base_dir=str(Path(cfg_path).parent),
        )

    def available_layers(self):
        layers = ['rotations']
        layers.extend(sorted(self._layer_files.keys()))
        return layers

    # Rotation
    def rotation_model(self) -> pygplates.RotationModel:
        if self._rotation_model is None:
            self._rotation_model = pygplates.RotationModel(self._rotation_file)
        return self._rotation_model

    # Layer paths from config['Layers']
    def feature(self, name: str) -> Union[str, List[str]]:
        """
        Return the resolved local filesystem path(s) for a named layer.

        If the JSON value is a single path (string) -> returns str.
        If the JSON value is a list of paths -> returns list[str].
        """
        if name not in self._layer_files:
            raise ValueError(f"Unknown layer '{name}'. Available: {sorted(self._layer_files.keys())}")
        spec = self._layer_files[name]
        if isinstance(spec, str):
            return _resolve_local(spec, self._base_dir)
        # spec is a list of string paths
        return [_resolve_local(p, self._base_dir) for p in spec]

    # Metadata accessors
    @property
    def big_time(self) -> Optional[float]:
        return self._big_time

    @property
    def small_time(self) -> Optional[float]:
        return self._small_time

    @property
    def time_dependent_rasters(self) -> Dict[str, object]:
        return self._time_dep_rasters

    # Convenience getters inspired by plate-model-manager API
    def _find_layer_key(self, candidate_names):
        # Case-insensitive and simple plural/singular handling
        lower_to_key = {k.lower(): k for k in self._layer_files.keys()}
        for name in candidate_names:
            key = lower_to_key.get(name.lower())
            if key:
                return key
        return None

    def get_rotation_file(self) -> str:
        """Return the resolved path to the rotations .rot file."""
        return self._rotation_file

    def get_coastlines(self) -> Union[str, List[str]]:
        key = self._find_layer_key(['Coastlines', 'Coastline'])
        if not key:
            raise ValueError("No 'Coastlines' layer defined in JSON 'Layers'.")
        return self.feature(key)

    def get_topologies(self) -> Union[str, List[str]]:
        key = self._find_layer_key(['Topologies', 'Topology'])
        if not key:
            raise ValueError("No 'Topologies' layer defined in JSON 'Layers'.")
        return self.feature(key)

    def get_COBs(self) -> Union[str, List[str]]:
        key = self._find_layer_key(['COBs', 'COB'])
        if not key:
            raise ValueError("No 'COBs' layer defined in JSON 'Layers'.")
        return self.feature(key)


def load_plate_model_from_json(config_path: str) -> ReconModel:
    '''Convenience function: build a ReconModel from a JSON config.'''
    return ReconModel.from_json(config_path)


# Helpers
def _is_url(path_str: str) -> bool:
    return path_str.startswith('http://') or path_str.startswith('https://')


def _local_path(path_str: str) -> str:
    if _is_url(path_str):
        raise ValueError(
            f"Remote URLs are not supported in this minimal loader: {path_str}. "
            "Please provide a local filesystem path."
        )
    return str(Path(path_str).expanduser())


def _resolve_local(path_str: str, base_dir: Optional[str]) -> str:
    if _is_url(path_str):
        raise ValueError(
            f"Remote URLs are not supported in this minimal loader: {path_str}. "
            "Please provide a local filesystem path."
        )
    p = Path(path_str).expanduser()
    if not p.is_absolute() and base_dir:
        p = Path(base_dir) / p
    return str(p.resolve())


def _assert_local_file(path_str: str, label: str, base_dir: Optional[str] = None) -> None:
    p = _resolve_local(path_str, base_dir)
    if not os.path.isfile(p):
        raise FileNotFoundError(f"{label} not found: {p}")


class iPlateModel(ReconModel):
    '''
    Convenience class so you can do:

        pm = iPlateModel("/abs/path/model.json")

    Equivalent to ReconModel.from_json(...).
    '''

    def __init__(self, config_json_path: str) -> None:
        # Parse JSON once, then call base constructor with resolved fields
        parsed = ReconModel.from_json(config_json_path)
        super().__init__(
            rotation_file=parsed._rotation_file,
            layers=parsed._layer_files,
            big_time=parsed._big_time,
            small_time=parsed._small_time,
            time_dep_rasters=parsed._time_dep_rasters,
            base_dir=parsed._base_dir,
        )


