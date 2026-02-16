#!/usr/bin/env python3
"""
OS OpenMap Local → DXF Converter
=================================
Converts OS OpenMap Local data (.gpkg or .gml) into AutoCAD-compatible .dxf.

Each layer becomes a DXF layer. Polygons → closed polylines (no hatches).
Lines → polylines. Points → point entities.

Usage:
    python gpkg2dxf.py <input.gpkg>              Single GeoPackage
    python gpkg2dxf.py <input.gml>               Single GML file
    python gpkg2dxf.py <folder_of_gml_files/>    Folder containing .gml files
    python gpkg2dxf.py                           Opens file picker (GUI)

Coordinates preserved as British National Grid (EPSG:27700) in metres.

Requirements:
    pip install geopandas ezdxf tqdm

DISCLAIMER
----------
This tool is provided "as is" without warranty of any kind, express or
implied. The output is derived from Ordnance Survey OpenData which is
subject to Crown copyright and the Open Government Licence. Users are
responsible for ensuring their use complies with the OS terms and
conditions: https://www.ordnancesurvey.co.uk/legal
The author accepts no liability for any consequences arising from use
of this tool or its output.

Author: Built for Jake White Architecture
Licence: MIT
"""

import sys
import os
import glob
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import geopandas as gpd
import ezdxf
import pyogrio
import numpy as np
from shapely import get_coordinates, get_parts
from shapely.geometry import (
    Polygon, MultiPolygon, LineString, MultiLineString,
    Point, MultiPoint, GeometryCollection,
)

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

try:
    import tkinter as tk
    from tkinter import filedialog
    HAS_TK = True
except ImportError:
    HAS_TK = False


# ── DXF layer colours (AutoCAD Color Index) ──────────────────────────────
LAYER_COLOURS = {
    "building": 252, "importantbuilding": 8,
    "road": 7, "roadtunnel": 9, "motorwayjunction": 5,
    "roundabout": 7, "railwaytrack": 8, "railwaytunnel": 8,
    "railwaystation": 8,
    "surfacewater_area": 4, "surfacewater_line": 4,
    "tidalwater": 4, "tidalboundary": 4, "foreshore": 42,
    "woodland": 3, "functionalsite": 62,
    "electricitytransmissionline": 1, "carchargingpoint": 2,
    "namedplace": 7, "roadname": 7,
}


def clean_layer_name(name):
    name = os.path.splitext(name)[0]
    parts = name.split("_", 1)
    if len(parts) == 2 and len(parts[0]) == 2 and parts[0].isupper():
        return parts[1]
    return name


def get_colour(layer_name):
    return LAYER_COLOURS.get(layer_name.lower(), 7)


# ─── FAST GEOMETRY → DXF (no iterrows) ───────────────────────────────────

def _polygon_to_rings(poly):
    """Extract all rings (exterior + interiors) as coordinate lists."""
    rings = []
    coords = poly.exterior.coords
    if len(coords) >= 3:
        rings.append([(c[0], c[1]) for c in coords])
    for interior in poly.interiors:
        coords = interior.coords
        if len(coords) >= 3:
            rings.append([(c[0], c[1]) for c in coords])
    return rings


def _geom_to_entities(geom):
    """
    Convert a single geometry to a list of (type, points, closed) tuples.
    Returns list of ('lwpoly', [(x,y),...], True/False) or ('point', (x,y), None).
    """
    if geom is None or geom.is_empty:
        return []

    results = []
    gt = geom.geom_type

    if gt == "Polygon":
        for ring in _polygon_to_rings(geom):
            results.append(("lwpoly", ring, True))

    elif gt == "MultiPolygon":
        for poly in geom.geoms:
            for ring in _polygon_to_rings(poly):
                results.append(("lwpoly", ring, True))

    elif gt == "LineString":
        coords = list(geom.coords)
        if len(coords) >= 2:
            results.append(("lwpoly", [(c[0], c[1]) for c in coords], False))

    elif gt == "MultiLineString":
        for line in geom.geoms:
            coords = list(line.coords)
            if len(coords) >= 2:
                results.append(
                    ("lwpoly", [(c[0], c[1]) for c in coords], False))

    elif gt == "Point":
        results.append(("point", (geom.x, geom.y), None))

    elif gt == "MultiPoint":
        for pt in geom.geoms:
            results.append(("point", (pt.x, pt.y), None))

    elif gt == "GeometryCollection":
        for sub in geom.geoms:
            results.extend(_geom_to_entities(sub))

    return results


def process_layer_to_entities(gdf, dxf_layer):
    """
    Convert an entire GeoDataFrame to a flat list of DXF-ready tuples.
    Runs geometry extraction in bulk — much faster than iterrows.
    Returns: list of (dxf_layer, type, points, closed)
    """
    entities = []
    geom_array = gdf.geometry.values  # shapely array

    for geom in geom_array:
        for etype, pts, closed in _geom_to_entities(geom):
            entities.append((dxf_layer, etype, pts, closed))

    return entities


# ─── CONCURRENT LAYER READING ────────────────────────────────────────────

def read_layer(source_path, internal_layer, dxf_layer):
    """
    Read a single layer from disk and convert to DXF entity tuples.
    Designed to run in a thread (GML/GPKG reading is I/O bound).
    Returns: (dxf_layer, feature_count, entity_list, error_or_None)
    """
    try:
        read_kwargs = {}
        if internal_layer:
            read_kwargs["layer"] = internal_layer
        gdf = gpd.read_file(source_path, **read_kwargs)
        feature_count = len(gdf)
        entities = process_layer_to_entities(gdf, dxf_layer)
        return (dxf_layer, feature_count, entities, None)
    except Exception as e:
        return (dxf_layer, 0, [], str(e))


# ─── SOURCE DISCOVERY ────────────────────────────────────────────────────

def discover_sources(input_path):
    """
    Return list of (file_path, internal_layer_name_or_None, dxf_layer) tuples.
    """
    sources = []

    if os.path.isdir(input_path):
        gml_files = sorted(set(
            glob.glob(os.path.join(input_path, "*.gml")) +
            glob.glob(os.path.join(input_path, "**/*.gml"), recursive=True)
        ))
        if not gml_files:
            print(f"ERROR: No .gml files found in {input_path}")
            sys.exit(1)
        for f in gml_files:
            sources.append((f, None, clean_layer_name(os.path.basename(f))))

    elif input_path.lower().endswith((".gpkg", ".gml")):
        try:
            layers = pyogrio.list_layers(input_path)
        except Exception as e:
            print(f"ERROR reading file: {e}")
            sys.exit(1)
        for name, geom_type in layers:
            if geom_type is None or geom_type == "":
                continue
            sources.append((input_path, name, clean_layer_name(name)))
    else:
        print(f"ERROR: Unrecognised file type: {input_path}")
        sys.exit(1)

    return sources


# ─── MAIN CONVERSION ─────────────────────────────────────────────────────

def convert_to_dxf(input_path, dxf_path=None, max_workers=4):
    """Main conversion with concurrent reading and batched DXF writing."""

    if not os.path.exists(input_path):
        print(f"ERROR: Path not found: {input_path}")
        sys.exit(1)

    if dxf_path is None:
        if os.path.isdir(input_path):
            dxf_path = os.path.join(
                input_path,
                os.path.basename(input_path.rstrip("/\\")) + ".dxf")
        else:
            dxf_path = os.path.splitext(input_path)[0] + ".dxf"

    if os.path.isfile(input_path):
        size_str = f"{os.path.getsize(input_path) / (1024*1024):.1f} MB"
    else:
        size_str = "folder"

    print(f"\n{'='*60}")
    print(f"  OS OpenMap Local → DXF Converter")
    print(f"{'='*60}")
    print(f"  Input:  {input_path} ({size_str})")
    print(f"  Output: {dxf_path}")
    print(f"  Workers: {max_workers} threads")
    print(f"{'='*60}")
    print(f"\n  DISCLAIMER: Output derived from OS OpenData.")
    print(f"  Subject to Crown copyright & Open Government Licence.")
    print(f"  See: ordnancesurvey.co.uk/legal")
    print(f"  No liability accepted. Use at your own risk.\n")

    # ── Discover sources ──────────────────────────────────────────────
    print("Scanning input data...")
    sources = discover_sources(input_path)
    print(f"Found {len(sources)} layers to process:\n")
    for _, _, dxf_layer in sources:
        print(f"  • {dxf_layer}")
    print()

    # ── Create DXF document ───────────────────────────────────────────
    doc = ezdxf.new("R2010")
    msp = doc.modelspace()

    # Pre-create all layers
    for _, _, dxf_layer in sources:
        if dxf_layer not in [l.dxf.name for l in doc.layers]:
            doc.layers.add(dxf_layer, color=get_colour(dxf_layer))

    # ── Read layers concurrently ──────────────────────────────────────
    start_time = time.time()
    all_entities = []
    total_features = 0
    completed = 0

    print(f"Reading layers ({max_workers} threads)...\n")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {}
        for source_path, internal_layer, dxf_layer in sources:
            f = executor.submit(
                read_layer, source_path, internal_layer, dxf_layer)
            futures[f] = (internal_layer or os.path.basename(source_path),
                          dxf_layer)

        for future in as_completed(futures):
            display_name, dxf_layer = futures[future]
            completed += 1
            layer_name, feat_count, entities, error = future.result()

            if error:
                print(f"  [{completed}/{len(sources)}] {display_name:<35} "
                      f"ERROR: {error}")
            else:
                all_entities.extend(entities)
                total_features += feat_count
                print(f"  [{completed}/{len(sources)}] {display_name:<35} "
                      f"{feat_count:>8,} features → "
                      f"{len(entities):>8,} entities")

    read_time = time.time() - start_time
    print(f"\n  Read complete: {total_features:,} features in {read_time:.1f}s")

    # ── Write entities to DXF (single-threaded, ezdxf isn't threadsafe)
    print(f"\n  Writing {len(all_entities):,} entities to DXF...")
    write_start = time.time()

    if HAS_TQDM and len(all_entities) > 1000:
        iterator = tqdm(
            all_entities,
            desc="  Writing DXF",
            unit="ent",
            bar_format="  {l_bar}{bar:35}{r_bar}",
        )
    else:
        iterator = all_entities

    for dxf_layer, etype, pts, closed in iterator:
        if etype == "lwpoly":
            msp.add_lwpolyline(pts, close=closed,
                               dxfattribs={"layer": dxf_layer})
        elif etype == "point":
            msp.add_point(pts, dxfattribs={"layer": dxf_layer})

    write_time = time.time() - write_start

    # ── Save ──────────────────────────────────────────────────────────
    print(f"\n  Saving file...")
    save_start = time.time()
    doc.saveas(dxf_path)
    save_time = time.time() - save_start

    total_time = time.time() - start_time
    dxf_size_mb = os.path.getsize(dxf_path) / (1024 * 1024)

    print(f"\n{'='*60}")
    print(f"  DONE!")
    print(f"  Total entities: {len(all_entities):,}")
    print(f"  Output size:    {dxf_size_mb:.1f} MB")
    print(f"  ──────────────────────────────────")
    print(f"  Read time:      {read_time:.1f}s")
    print(f"  Write time:     {write_time:.1f}s")
    print(f"  Save time:      {save_time:.1f}s")
    print(f"  Total time:     {total_time:.1f}s")
    print(f"{'='*60}")
    print(f"\n  Open {os.path.basename(dxf_path)} in BricsCAD.")
    print(f"  Coordinates are British National Grid (EPSG:27700).")
    print(f"  Type ZOOM E to see everything.\n")


# ─── FILE PICKER ──────────────────────────────────────────────────────────

def pick_files():
    root = tk.Tk()
    root.withdraw()

    input_path = filedialog.askopenfilename(
        title="Select OS OpenMap Local file (.gpkg or .gml)",
        filetypes=[
            ("Geo files", "*.gpkg *.gml"),
            ("GeoPackage", "*.gpkg"),
            ("GML files", "*.gml"),
            ("All files", "*.*"),
        ],
    )
    if not input_path:
        input_path = filedialog.askdirectory(
            title="Or select a FOLDER containing .gml files")
        if not input_path:
            print("No input selected. Exiting.")
            sys.exit(0)

    if os.path.isdir(input_path):
        default_out = os.path.join(
            input_path,
            os.path.basename(input_path.rstrip("/\\")) + ".dxf")
    else:
        default_out = os.path.splitext(input_path)[0] + ".dxf"

    output_path = filedialog.asksaveasfilename(
        title="Save DXF As",
        initialfile=os.path.basename(default_out),
        initialdir=os.path.dirname(default_out),
        defaultextension=".dxf",
        filetypes=[("DXF files", "*.dxf"), ("All files", "*.*")],
    )
    if not output_path:
        output_path = default_out

    root.destroy()
    return input_path, output_path


# ─── ENTRY POINT ──────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Parse args
    workers = 8
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    for a in sys.argv[1:]:
        if a.startswith("--workers="):
            workers = int(a.split("=")[1])

    if args:
        input_path = args[0]
        output_path = args[1] if len(args) > 1 else None
    else:
        if not HAS_TK:
            print("Usage: python gpkg2dxf.py <input> [output.dxf] [--workers=4]")
            print("\n  <input> can be:")
            print("    - A .gpkg file (GeoPackage)")
            print("    - A .gml file (multiple internal layers)")
            print("    - A folder containing .gml files")
            print("\n  --workers=N  Number of threads (default 8)")
            print("\nOr just double-click for a file picker.")
            sys.exit(1)
        input_path, output_path = pick_files()

    convert_to_dxf(input_path, output_path, max_workers=workers)
