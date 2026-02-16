# os2dxf

Convert Ordnance Survey OpenData into AutoCAD-compatible DXF files.

Built for architects and technologists who need OS mapping data in their CAD environment without paying for expensive conversion tools.

## What it does

- Reads **OS OpenMap Local** data (`.gml` or `.gpkg` format)
- Outputs a clean `.dxf` with each OS layer (Building, Road, Woodland, etc.) on its own DXF layer
- Polygons become closed polylines (outlines only, no hatches)
- Lines become polylines
- Points become point entities
- Coordinates preserved as **British National Grid (EPSG:27700)** in metres
- Multi-threaded reading for large files (default 8 threads)
- Progress bar via tqdm

## Quick start

```bash
pip install geopandas ezdxf tqdm
python os2dxf.py
```

Double-click or run without arguments to get a file picker. Or from the command line:

```bash
python os2dxf.py TQ.gml                       # single GML (auto-detects internal layers)
python os2dxf.py TQ.gml output.dxf             # specify output path
python os2dxf.py openmap_local.gpkg             # GeoPackage input
python os2dxf.py /path/to/gml_folder/           # folder of .gml files
python os2dxf.py TQ.gml --workers=16            # more threads
```

## Getting the data

1. Go to [OS Data Hub](https://osdatahub.os.uk/downloads/open/OpenMapLocal)
2. Select **OS OpenMap Local**
3. Choose your grid square and **GML** or **GeoPackage** format
4. Download and unzip
5. Point `os2dxf.py` at the `.gml` or `.gpkg` file

## DXF layer mapping

| OS Layer | DXF Layer | ACI Colour |
|---|---|---|
| Building | Building | 252 (light grey) |
| ImportantBuilding | ImportantBuilding | 8 (dark grey) |
| Road | Road | 7 (white) |
| RailwayTrack | RailwayTrack | 8 (dark grey) |
| SurfaceWater_Area | SurfaceWater_Area | 4 (cyan) |
| SurfaceWater_Line | SurfaceWater_Line | 4 (cyan) |
| Woodland | Woodland | 3 (green) |
| FunctionalSite | FunctionalSite | 62 (brown) |
| Foreshore | Foreshore | 42 (tan) |
| TidalWater | TidalWater | 4 (cyan) |
| ElectricityTransmissionLine | ElectricityTransmissionLine | 1 (red) |

Other layers default to colour 7 (white/black depending on background).

## In your CAD software

Open the `.dxf` and type `ZOOM E` to extents. Tested with BricsCAD — should work with any DXF-compatible CAD (AutoCAD, LibreCAD, FreeCAD, etc.).

## Requirements

- Python 3.10+
- geopandas
- ezdxf
- tqdm (optional, for progress bars)
- tkinter (included with most Python installs, for file picker)

## Disclaimer

This tool is provided "as is" without warranty of any kind. Output is derived from Ordnance Survey OpenData, subject to Crown copyright and the [Open Government Licence](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/). Users are responsible for ensuring their use complies with [OS terms and conditions](https://www.ordnancesurvey.co.uk/legal). The author accepts no liability for any consequences arising from use of this tool or its output.

## Licence

MIT
