"""
sinkhole_detect.py
------------------
LiDAR-based sinkhole detection from a DEM GeoTIFF.

Compatible with both Jupyter Notebook and command-line execution.

JUPYTER USAGE:
    - Edit the JUPYTER_CONFIG block at the top of this file
    - Run the script as a cell or with %run sinkhole_detect.py

COMMAND LINE USAGE:
    python sinkhole_detect.py --input Lidarsample.tif --output results/
    python sinkhole_detect.py --input Lidarsample.tif --min-depth 0.5 --sigma 1.0

OUTPUTS (written to output directory):
    - *_sinkholes.tif      RGB composite raster with highlighted sinkholes
    - *_sinkholes.geojson  Vector point file with sinkhole attributes
    - *_sinkholes.csv      CSV log of all detected sinkholes
    - *_preview.png        Visualization image with legend

DEPENDENCIES:
    pip install richdem rasterio numpy scipy scikit-image matplotlib geopandas shapely
"""

import csv
import os
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')


# ═══════════════════════════════════════════════════════════════
#  JUPYTER / DIRECT CONFIGURATION
#  If running in Jupyter, edit these values and run the script.
#  These are ignored when running from the command line.
# ═══════════════════════════════════════════════════════════════

JUPYTER_CONFIG = {
    'input_file':                   ,           # Path to your DEM GeoTIFF
    'output_dir':         'results/',           # Output folder (created if missing)
    'min_sinkhole_depth': 15,                  # Minimum depression depth to detect
    'max_sinkhole_depth': 100.0,                 # Maximum depression depth (filters artifacts)
    'min_area':           10,                   # Minimum region size in pixels
    'min_circularity':    0.4,                  # Circularity threshold 0-1 (1 = perfect circle)
    'max_axis_ratio':     4.5,                  # Filters elongated/linear features like ditches
    'gaussian_sigma':     0.7,                  # Smoothing strength (0 = disabled)
}

# Gradient colors: shallow edge = yellow-orange, deep center = black
SHALLOW_COLOR = np.array([255, 165, 0], dtype=np.float32)
DEEP_COLOR    = np.array([0,   0,   0], dtype=np.float32)


# ═══════════════════════════════════════════════════════════════
#  STEP 1: LOAD DEM
# ═══════════════════════════════════════════════════════════════

def load_dem(filepath: str):
    """
    Load a single-band elevation GeoTIFF.
    Returns: dem (float32 ndarray), profile, affine transform, CRS.
    Raises FileNotFoundError or ValueError on bad input.
    """
    import rasterio

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Input file not found: {filepath}")

    with rasterio.open(filepath) as src:
        dem_ma    = src.read(1, masked=True)
        profile   = src.profile.copy()
        transform = src.transform
        crs       = src.crs

    dem = dem_ma.astype(np.float32).data
    dem[dem_ma.mask] = np.nan
    dem[dem < 0]     = np.nan  # Remove invalid negative elevations

    if np.all(np.isnan(dem)):
        raise ValueError(
            "DEM contains only NaN values. "
            "Check that your input is a valid single-band elevation raster."
        )

    print(f"  Loaded : {filepath}")
    print(f"  Shape  : {dem.shape[0]} rows x {dem.shape[1]} cols")
    print(f"  CRS    : {crs}")
    print(f"  Elev   : {np.nanmin(dem):.2f} to {np.nanmax(dem):.2f}")

    return dem, profile, transform, crs


# ═══════════════════════════════════════════════════════════════
#  STEP 2: DETECT DEPRESSIONS
# ═══════════════════════════════════════════════════════════════

def detect_depressions(dem: np.ndarray, config: dict) -> np.ndarray:
    """
    Optional Gaussian smoothing -> RichDEM depression fill -> depth difference.
    Returns a depression-depth array with values outside min/max set to 0.
    """
    import richdem as rd
    from scipy import ndimage as ndi

    sigma = config['gaussian_sigma']
    if sigma > 0:
        print(f"  Applying Gaussian smoothing (sigma={sigma})...")
        smoothed = ndi.gaussian_filter(np.where(np.isnan(dem), 0, dem), sigma=sigma)
        smoothed[np.isnan(dem)] = np.nan
    else:
        smoothed = dem.copy()

    print("  Filling depressions with RichDEM...")
    rdem   = rd.rdarray(smoothed, no_data=np.nan)
    filled = rd.FillDepressions(rdem, in_place=False)
    diff   = np.array(filled) - np.array(rdem)

    depression = diff.copy()
    depression[depression < config['min_sinkhole_depth']] = 0
    depression[depression > config['max_sinkhole_depth']] = 0

    print(f"  Depression pixels found: {np.count_nonzero(depression):,}")
    return depression


# ═══════════════════════════════════════════════════════════════
#  STEP 3: FILTER REGIONS
# ═══════════════════════════════════════════════════════════════

def filter_regions(depression: np.ndarray, config: dict):
    """
    Label connected depression regions and filter by:
      - Minimum area in pixels
      - Minimum circularity  (rejects irregular noise)
      - Maximum axis ratio   (rejects ditches and linear features)

    Returns: (labeled array, list of (region, circularity, axis_ratio) tuples)
    """
    from scipy import ndimage as ndi
    from skimage import measure

    binary     = (depression > 0).astype(np.uint8)
    labeled, _ = ndi.label(binary)
    regions    = measure.regionprops(labeled, intensity_image=depression)

    qualified = []
    rejected  = {'area': 0, 'circularity': 0, 'axis_ratio': 0}

    for region in regions:
        # Area filter
        if region.area < config['min_area']:
            rejected['area'] += 1
            continue

        # Circularity filter
        if region.perimeter == 0:
            continue
        circularity = 4 * np.pi * region.area / (region.perimeter ** 2)
        if circularity < config['min_circularity']:
            rejected['circularity'] += 1
            continue

        # Axis ratio filter - rejects elongated linear features
        if region.minor_axis_length > 0:
            axis_ratio = region.major_axis_length / region.minor_axis_length
        else:
            axis_ratio = float('inf')
        if axis_ratio > config['max_axis_ratio']:
            rejected['axis_ratio'] += 1
            continue

        qualified.append((region, circularity, axis_ratio))
        print(f"    + Region {region.label:>4}  "
              f"area={region.area:>6}px  "
              f"circ={circularity:.2f}  "
              f"axis={axis_ratio:.2f}  "
              f"max_depth={region.max_intensity:.3f}")

    print(f"\n  Rejected - too small : {rejected['area']}")
    print(f"  Rejected - low circ  : {rejected['circularity']}")
    print(f"  Rejected - elongated : {rejected['axis_ratio']}")
    print(f"  Qualified sinkholes  : {len(qualified)}")

    return labeled, qualified


# ═══════════════════════════════════════════════════════════════
#  STEP 4: BUILD RGB OVERLAY  (fully vectorized)
# ═══════════════════════════════════════════════════════════════

def build_overlay(qualified_regions, depression: np.ndarray,
                  labeled: np.ndarray, dem: np.ndarray) -> np.ndarray:
    """
    Build a 3-channel RGB composite using vectorized NumPy operations.
    Sinkhole pixels: shallow edge -> yellow-orange, deep center -> black.
    All other pixels: grayscale DEM background.
    """
    overlay = np.zeros((3, dem.shape[0], dem.shape[1]), dtype=np.uint8)

    for region, _, _ in qualified_regions:
        region_mask  = (labeled == region.label)
        region_depth = depression[region_mask]

        depth_range = region_depth.ptp()
        if depth_range < 1e-6:
            normalized = np.zeros_like(region_depth)
        else:
            normalized = (region_depth - region_depth.min()) / (depth_range + 1e-6)
        normalized = np.clip(normalized, 0, 1)

        # Vectorized interpolation: shallow=yellow-orange, deep=black
        color_array = (SHALLOW_COLOR[None, :] * (1 - normalized[:, None])).astype(np.uint8)

        coords = np.where(region_mask)
        overlay[0][coords] = color_array[:, 0]
        overlay[1][coords] = color_array[:, 1]
        overlay[2][coords] = color_array[:, 2]

    # Grayscale DEM background
    gray       = np.interp(dem, (np.nanmin(dem), np.nanmax(dem)), (0, 255)).astype(np.uint8)
    background = np.stack([gray] * 3, axis=0)
    final      = np.where(overlay.sum(axis=0, keepdims=True) > 0, overlay, background)

    return final


# ═══════════════════════════════════════════════════════════════
#  STEP 5: SAVE RASTER OUTPUT
# ═══════════════════════════════════════════════════════════════

def save_raster(final: np.ndarray, profile: dict, output_path: str):
    """Save the 3-band RGB composite as a GeoTIFF."""
    import rasterio
    profile.update({'count': 3, 'dtype': 'uint8'})
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(final)
    print(f"  Raster  -> {output_path}")


# ═══════════════════════════════════════════════════════════════
#  STEP 6: SAVE VECTOR + CSV OUTPUTS
# ═══════════════════════════════════════════════════════════════

def save_vector_outputs(qualified_regions, depression: np.ndarray,
                        labeled: np.ndarray, transform, crs,
                        geojson_path: str, csv_path: str):
    """
    Write sinkhole centroids to GeoJSON and CSV.
    Attributes: id, area_px, max_depth, mean_depth, circularity,
                axis_ratio, centroid_x, centroid_y (real-world CRS).
    """
    import geopandas as gpd
    from rasterio.transform import xy
    from shapely.geometry import Point

    if not qualified_regions:
        print("  No sinkholes detected - vector outputs skipped.")
        return

    records = []
    for i, (region, circularity, axis_ratio) in enumerate(qualified_regions, start=1):
        row, col = region.centroid
        cx, cy   = xy(transform, row, col)
        records.append({
            'id':          i,
            'area_px':     region.area,
            'max_depth':   round(float(region.max_intensity), 3),
            'mean_depth':  round(float(region.mean_intensity), 3),
            'circularity': round(circularity, 3),
            'axis_ratio':  round(axis_ratio, 3),
            'centroid_x':  round(cx, 4),
            'centroid_y':  round(cy, 4),
            'geometry':    Point(cx, cy),
        })

    # GeoJSON
    gdf = gpd.GeoDataFrame(records, geometry='geometry', crs=crs)
    gdf.to_file(geojson_path, driver='GeoJSON')
    print(f"  GeoJSON -> {geojson_path}")

    # CSV
    csv_records = [{k: v for k, v in r.items() if k != 'geometry'} for r in records]
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=csv_records[0].keys())
        writer.writeheader()
        writer.writerows(csv_records)
    print(f"  CSV     -> {csv_path}")


# ═══════════════════════════════════════════════════════════════
#  STEP 7: VISUALIZE
# ═══════════════════════════════════════════════════════════════

def visualize(final: np.ndarray, output_path: str):
    """Display and save the composite with a color legend."""
    from matplotlib.patches import Patch

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(np.moveaxis(final, 0, -1))
    ax.set_title("Detected Sinkholes - Circularity & Axis-Ratio Filtered", fontsize=13)
    ax.axis("off")

    legend_elements = [
        Patch(facecolor='#FFA500', label='Shallow depression edge'),
        Patch(facecolor='#1a1a1a', label='Deep depression center'),
        Patch(facecolor='#808080', label='DEM background (grayscale)'),
    ]
    ax.legend(handles=legend_elements, loc='lower right',
              framealpha=0.85, fontsize=10)

    plt.tight_layout()
    preview_path = output_path.replace('.tif', '_preview.png')
    plt.savefig(preview_path, dpi=150, bbox_inches='tight')
    print(f"  Preview -> {preview_path}")
    plt.show()


# ═══════════════════════════════════════════════════════════════
#  CORE RUNNER  (shared by both Jupyter and CLI)
# ═══════════════════════════════════════════════════════════════

def run(input_file: str, output_dir: str, config: dict):
    """
    Execute the full sinkhole detection pipeline.
    Called by both Jupyter and CLI entry points.
    """
    os.makedirs(output_dir, exist_ok=True)
    base        = os.path.splitext(os.path.basename(input_file))[0]
    raster_out  = os.path.join(output_dir, f"{base}_sinkholes.tif")
    geojson_out = os.path.join(output_dir, f"{base}_sinkholes.geojson")
    csv_out     = os.path.join(output_dir, f"{base}_sinkholes.csv")

    divider = '-' * 55

    print(f"\n{divider}")
    print(f"  Sinkhole Detection  |  {base}")
    print(divider)

    try:
        print("\n[1/7] Loading DEM...")
        dem, profile, transform, crs = load_dem(input_file)

        print("\n[2/7] Detecting depressions...")
        depression = detect_depressions(dem, config)

        print("\n[3/7] Filtering regions...")
        labeled, qualified = filter_regions(depression, config)

        print("\n[4/7] Building RGB overlay...")
        final = build_overlay(qualified, depression, labeled, dem)

        print("\n[5/7] Saving raster...")
        save_raster(final, profile, raster_out)

        print("\n[6/7] Saving vector outputs...")
        save_vector_outputs(qualified, depression, labeled,
                            transform, crs, geojson_out, csv_out)

        print("\n[7/7] Visualizing...")
        visualize(final, raster_out)

        print(f"\n{divider}")
        print(f"  Done. {len(qualified)} sinkhole(s) detected.")
        print(f"{divider}\n")

    except (FileNotFoundError, ValueError) as e:
        print(f"\n  ERROR: {e}\n", file=sys.stderr)
        sys.exit(1)


# ═══════════════════════════════════════════════════════════════
#  ENVIRONMENT DETECTION  - Jupyter vs. Command Line
# ═══════════════════════════════════════════════════════════════

def _is_jupyter() -> bool:
    """Return True if running inside a Jupyter / IPython kernel."""
    try:
        shell = get_ipython().__class__.__name__
        return shell in ('ZMQInteractiveShell', 'TerminalInteractiveShell')
    except NameError:
        return False


def _parse_args():
    """Parse command-line arguments (only used outside Jupyter)."""
    import argparse
    parser = argparse.ArgumentParser(
        description='Detect sinkholes in a LiDAR DEM GeoTIFF.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input',     required=True,
                        help='Path to input DEM GeoTIFF')
    parser.add_argument('--output',    default='results/',
                        help='Output directory')
    parser.add_argument('--min-depth', type=float,
                        default=JUPYTER_CONFIG['min_sinkhole_depth'],
                        help='Minimum depression depth')
    parser.add_argument('--max-depth', type=float,
                        default=JUPYTER_CONFIG['max_sinkhole_depth'],
                        help='Maximum depression depth')
    parser.add_argument('--min-area',  type=int,
                        default=JUPYTER_CONFIG['min_area'],
                        help='Minimum region area in pixels')
    parser.add_argument('--min-circ',  type=float,
                        default=JUPYTER_CONFIG['min_circularity'],
                        help='Minimum circularity (0-1)')
    parser.add_argument('--max-axis',  type=float,
                        default=JUPYTER_CONFIG['max_axis_ratio'],
                        help='Maximum axis ratio (filters linear features)')
    parser.add_argument('--sigma',     type=float,
                        default=JUPYTER_CONFIG['gaussian_sigma'],
                        help='Gaussian smoothing sigma (0 = disabled)')
    return parser.parse_args()


# ═══════════════════════════════════════════════════════════════
#  ENTRY POINT
# ═══════════════════════════════════════════════════════════════

if _is_jupyter():
    # Jupyter: use JUPYTER_CONFIG defined at the top of this file
    print("Running in Jupyter mode - using JUPYTER_CONFIG settings.")
    run(
        input_file = JUPYTER_CONFIG['input_file'],
        output_dir = JUPYTER_CONFIG['output_dir'],
        config     = {k: v for k, v in JUPYTER_CONFIG.items()
                      if k not in ('input_file', 'output_dir')}
    )

else:
    # Command line: parse arguments
    if __name__ == '__main__':
        args = _parse_args()
        run(
            input_file = args.input,
            output_dir = args.output,
            config     = {
                'min_sinkhole_depth': args.min_depth,
                'max_sinkhole_depth': args.max_depth,
                'min_area':           args.min_area,
                'min_circularity':    args.min_circ,
                'max_axis_ratio':     args.max_axis,
                'gaussian_sigma':     args.sigma,
            }
        )
