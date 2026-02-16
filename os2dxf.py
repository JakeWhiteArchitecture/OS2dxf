#!/usr/bin/env python3
"""
os2dxf — OS OpenMap Local → DXF Converter
==========================================
Converts OS OpenMap Local data (.gml or .gpkg) into per-layer DXF files.

Each OS layer (Building, Road, Woodland, etc.) becomes its own .dxf file
in an output folder. Fully parallel: each thread reads a layer from the
source AND writes its own DXF — no shared file bottleneck.

Use XREF in BricsCAD/AutoCAD to attach the layers you need.

Usage:
    python os2dxf.py                             File picker (GUI)
    python os2dxf.py TQ.gml                      Single GML
    python os2dxf.py data.gpkg                    GeoPackage
    python os2dxf.py /folder/of/gmls/             Folder of .gml files
    python os2dxf.py TQ.gml --workers=16          More threads

Output:
    Creates a folder (e.g. TQ_dxf/) containing one .dxf per layer:
        Building.dxf, Road.dxf, Woodland.dxf, etc.

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


# ── Layer colours (AutoCAD Color Index) ───────────────────────────────────
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


# ── Geometry → DXF entities ───────────────────────────────────────────────

def _polygon_to_rings(poly):
    rings = []
    coords = poly.exterior.coords
    if len(coords) >= 3:
        rings.append([(c[0], c[1]) for c in coords])
    for interior in poly.interiors:
        coords = interior.coords
        if len(coords) >= 3:
            rings.append([(c[0], c[1]) for c in coords])
    return rings


def add_geometry(msp, geom, layer_name):
    """Add a geometry to the modelspace. Returns entity count."""
    if geom is None or geom.is_empty:
        return 0

    gt = geom.geom_type
    count = 0
    attribs = {"layer": layer_name}

    if gt == "Polygon":
        for ring in _polygon_to_rings(geom):
            msp.add_lwpolyline(ring, close=True, dxfattribs=attribs)
            count += 1

    elif gt == "MultiPolygon":
        for poly in geom.geoms:
            for ring in _polygon_to_rings(poly):
                msp.add_lwpolyline(ring, close=True, dxfattribs=attribs)
                count += 1

    elif gt == "LineString":
        coords = list(geom.coords)
        if len(coords) >= 2:
            msp.add_lwpolyline(
                [(c[0], c[1]) for c in coords], close=False,
                dxfattribs=attribs)
            count += 1

    elif gt == "MultiLineString":
        for line in geom.geoms:
            coords = list(line.coords)
            if len(coords) >= 2:
                msp.add_lwpolyline(
                    [(c[0], c[1]) for c in coords], close=False,
                    dxfattribs=attribs)
                count += 1

    elif gt == "Point":
        msp.add_point((geom.x, geom.y), dxfattribs=attribs)
        count += 1

    elif gt == "MultiPoint":
        for pt in geom.geoms:
            msp.add_point((pt.x, pt.y), dxfattribs=attribs)
            count += 1

    elif gt == "GeometryCollection":
        for sub in geom.geoms:
            count += add_geometry(msp, sub, layer_name)

    return count


# ── Per-layer worker (read + write in one thread) ─────────────────────────

def process_layer(source_path, internal_layer, dxf_layer, output_dir):
    """
    Complete pipeline for one layer: read from source → write .dxf.
    Each thread runs this independently — no shared state.
    Returns: (dxf_layer, feature_count, entity_count, dxf_path, elapsed, error)
    """
    t0 = time.time()
    dxf_path = os.path.join(output_dir, f"{dxf_layer}.dxf")

    try:
        # Read
        read_kwargs = {}
        if internal_layer:
            read_kwargs["layer"] = internal_layer
        gdf = gpd.read_file(source_path, **read_kwargs)
        feature_count = len(gdf)

        if feature_count == 0:
            return (dxf_layer, 0, 0, None, time.time() - t0, None)

        # Create DXF
        doc = ezdxf.new("R2010")
        colour = get_colour(dxf_layer)
        doc.layers.add(dxf_layer, color=colour)
        msp = doc.modelspace()

        # Write entities
        entity_count = 0
        for geom in gdf.geometry.values:
            entity_count += add_geometry(msp, geom, dxf_layer)

        # Save
        doc.saveas(dxf_path)

        return (dxf_layer, feature_count, entity_count, dxf_path,
                time.time() - t0, None)

    except Exception as e:
        return (dxf_layer, 0, 0, None, time.time() - t0, str(e))


# ── Source discovery ──────────────────────────────────────────────────────

def discover_sources(input_path):
    """Return list of (file_path, internal_layer_or_None, dxf_layer) tuples."""
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


# ── Main ──────────────────────────────────────────────────────────────────

def convert_to_dxf(input_path, output_dir=None, max_workers=8):
    if not os.path.exists(input_path):
        print(f"ERROR: Path not found: {input_path}")
        sys.exit(1)

    # Output folder
    if output_dir is None:
        if os.path.isdir(input_path):
            base = os.path.basename(input_path.rstrip("/\\"))
        else:
            base = os.path.splitext(os.path.basename(input_path))[0]
        output_dir = os.path.join(os.path.dirname(os.path.abspath(input_path)),
                                  f"{base}_dxf")

    os.makedirs(output_dir, exist_ok=True)

    if os.path.isfile(input_path):
        size_str = f"{os.path.getsize(input_path) / (1024*1024):.0f} MB"
    else:
        size_str = "folder"

    print(f"\n{'='*60}")
    print(f"  os2dxf — OS OpenMap Local → DXF")
    print(f"{'='*60}")
    print(f"  Input:   {input_path} ({size_str})")
    print(f"  Output:  {output_dir}/")
    print(f"  Workers: {max_workers} threads")
    print(f"{'='*60}")
    print(f"\n  DISCLAIMER: Output derived from OS OpenData.")
    print(f"  Subject to Crown copyright & Open Government Licence.")
    print(f"  See: ordnancesurvey.co.uk/legal\n")

    # Discover layers
    print("  Scanning input data...")
    sources = discover_sources(input_path)
    print(f"  Found {len(sources)} layers:\n")
    for _, _, dxf_layer in sources:
        print(f"    • {dxf_layer}")

    # Process all layers in parallel
    print(f"\n  Processing ({max_workers} threads)...\n")
    start_time = time.time()

    total_features = 0
    total_entities = 0
    results = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {}
        for source_path, internal_layer, dxf_layer in sources:
            f = executor.submit(
                process_layer, source_path, internal_layer,
                dxf_layer, output_dir)
            futures[f] = dxf_layer

        if HAS_TQDM:
            pbar = tqdm(
                total=len(sources),
                desc="  Layers",
                unit="layer",
                bar_format="  {l_bar}{bar:30}{r_bar}",
            )

        for future in as_completed(futures):
            layer, feats, ents, path, elapsed, error = future.result()

            if error:
                print(f"    ✗ {layer:<30} ERROR: {error}")
            elif path:
                dxf_mb = os.path.getsize(path) / (1024 * 1024)
                print(f"    ✓ {layer:<30} "
                      f"{feats:>8,} feat → {ents:>8,} ent  "
                      f"({dxf_mb:>6.1f} MB, {elapsed:>5.1f}s)")
                total_features += feats
                total_entities += ents
                results.append((layer, path, dxf_mb))
            else:
                print(f"    – {layer:<30} empty, skipped")

            if HAS_TQDM:
                pbar.update(1)

        if HAS_TQDM:
            pbar.close()

    total_time = time.time() - start_time
    total_dxf_mb = sum(r[2] for r in results)

    print(f"\n{'='*60}")
    print(f"  DONE!")
    print(f"  Layers:         {len(results)}")
    print(f"  Total features: {total_features:,}")
    print(f"  Total entities: {total_entities:,}")
    print(f"  Total DXF size: {total_dxf_mb:.1f} MB")
    print(f"  Total time:     {total_time:.1f}s")
    print(f"  Output folder:  {output_dir}")
    print(f"{'='*60}")
    print(f"\n  In BricsCAD:")
    print(f"    1. Open any DXF as your base")
    print(f"    2. XREF → Attach the others")
    print(f"    3. ZOOM E")
    print(f"  Coordinates: British National Grid (EPSG:27700)\n")

    # Write a handy file list
    manifest = os.path.join(output_dir, "_layers.txt")
    with open(manifest, "w") as f:
        f.write(f"os2dxf output — {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"Source: {input_path}\n\n")
        for layer, path, mb in sorted(results):
            f.write(f"{os.path.basename(path):<40} {mb:>8.1f} MB\n")
        f.write(f"\nTotal: {total_dxf_mb:.1f} MB across {len(results)} files\n")


# ── File picker ───────────────────────────────────────────────────────────

def pick_input():
    root = tk.Tk()
    root.withdraw()

    input_path = filedialog.askopenfilename(
        title="Select OS OpenMap Local file (.gml or .gpkg)",
        filetypes=[
            ("Geo files", "*.gpkg *.gml"),
            ("GML files", "*.gml"),
            ("GeoPackage", "*.gpkg"),
            ("All files", "*.*"),
        ],
    )
    if not input_path:
        input_path = filedialog.askdirectory(
            title="Or select a FOLDER containing .gml files")
        if not input_path:
            print("No input selected. Exiting.")
            sys.exit(0)

    root.destroy()
    return input_path


# ── Entry point ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    workers = 8
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    for a in sys.argv[1:]:
        if a.startswith("--workers="):
            workers = int(a.split("=")[1])

    if args:
        input_path = args[0]
        output_dir = args[1] if len(args) > 1 else None
    else:
        if not HAS_TK:
            print("Usage: python os2dxf.py <input> [output_folder] [--workers=8]")
            print("\n  <input> can be:")
            print("    - A .gml file (single file with internal layers)")
            print("    - A .gpkg file (GeoPackage)")
            print("    - A folder containing .gml files")
            print("\n  Output: creates <name>_dxf/ folder with one .dxf per layer")
            print("\nOr double-click for a file picker.")
            sys.exit(1)
        input_path = pick_input()
        output_dir = None

    convert_to_dxf(input_path, output_dir, max_workers=workers)
