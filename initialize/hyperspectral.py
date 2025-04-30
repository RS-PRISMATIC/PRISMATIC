import logging
import os
import sys
# import whitebox
import zipfile
import ray

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Patch
import rasterio
from rasterio.plot import show
from rasterio.features import geometry_mask
import fiona
from sklearn.preprocessing import OneHotEncoder
# from keras.models import Sequential
# from keras.layers import Conv2D, MaxPooling2D, Flatten, Dense
import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry import mapping
from shapely.geometry import box
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from collections import Counter
import seaborn as sns

# data preparation
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, accuracy_score, precision_score, recall_score, f1_score, classification_report

# model selection
from sklearn.model_selection import (
    GridSearchCV,
    GroupKFold,
    GroupShuffleSplit,
    KFold,
    ShuffleSplit,
    StratifiedGroupKFold,
    StratifiedKFold,
    StratifiedShuffleSplit,
    TimeSeriesSplit,
    train_test_split,
)

# model evaluation
from sklearn.metrics import f1_score

# machine learning
from sklearn.ensemble import RandomForestClassifier as RF

# Save files
import joblib
from joblib import dump, load

from utils.download_functions import download_aop_files
from utils.apply_brdf_corrections import move_red_yellow_subset, \
                                        convert_hdf_to_envi, \
                                        create_config, \
                                        implement_brdf_correction , \
                                        convert_envi_to_tif
from pathlib import Path

# add environ
conda_env_path = Path(sys.executable).parent.parent
os.environ['PROJ_LIB'] = str(conda_env_path/'share'/'proj')

log = logging.getLogger(__name__)

def download_hyperspectral(site, year, data_raw_aop_path, hs_type):
    """Download raw, uncorrected hyperspectral data (raw and raster) for site-year
    Analog to 04-download_aop_imagery.R from https://github.com/earthlab/neon-veg

    Parameters
    ----------
    site : str
        Site name
    year : str
        Year of aop collection
    data_raw_aop_path : str
        Path to store the downloaded data
    hs_type : str
        "tile" (L3 1km x 1km tile) OR "flightline" (L1 flightline)

    Returns
    -------
    (str, str)
        Path to the result folder for hyperspectral uncorrected imagery (L3 tiles) ('*/hs_tile' or '*/hs_flightline')
        and the result raster folder ('*/tif')
    """
    path = Path(data_raw_aop_path)/site/year

    product_code = 'DP3.30010.001' #rgb
    file_type = 'tif' 
    p = path/file_type
    download_aop_files(product_code,
                        site,
                        year,
                        str(p),
                        match_string=file_type,
                        check_size=False)

    product_code = 'DP3.30026.001' #veg indices
    download_aop_files(product_code,
                       site,
                       year,
                       str(p),
                       match_string='.zip',
                       check_size=False)
    zip_files = [file for file in os.listdir(p) if file.endswith('.zip')] # get the list of files
    for zip_file in zip_files:  #for each zipfile
        with zipfile.ZipFile(Path(p/zip_file)) as item: # treat the file as a zip
            item.extractall(p)  # extract it in the working directory
            item.close()
    vi_error_files = [file for file in os.listdir(p) if file.endswith('_error.tif')]
    for vi_error_file in vi_error_files: 
        os.remove(Path(p/vi_error_file)) # remove error files

    if hs_type=="tile":
        product_code = 'DP3.30006.001'  #DP3.30006.002
        p = path/'hs_tile'
    elif hs_type=="flightline": 
        product_code = 'DP1.30006.001' #DP1.30006.002
        p = path/'hs_flightline'
    else:
        print("must specify hs_type argument")
    file_type = 'h5'
    download_aop_files(product_code,
                       site,
                       year,
                       str(p),
                       match_string=file_type,
                       check_size=False)
    
    return str(p)

# def generate_pft_reference(sites, data_raw_inv_path, data_int_path, trait_table_path):
#     """ Generate a csv that translates each species into a PFT usable as training data
#     """
#     r_source = ro.r['source']
#     r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

#     # Generate PFT reference
#     generate_pft_reference = ro.r('generate_pft_reference')
#     Path(data_int_path).mkdir(parents=True, exist_ok=True)
#     pft_reference_path = generate_pft_reference(sites, data_raw_inv_path, data_int_path, trait_table_path)
#     log.info('Generated PFT reference for: {sites}'
#              f'{pft_reference_path}')
#     #ais now where to feed in pft_reference csv path


def correct_flightlines(site, year_inv, year_aop, data_raw_aop_path, data_int_path):
    """Correct raw L1 NEON AOP flightlines.
        1) Apply BRDF/topographic corrections
            output format: *.envi 
        2)  Convert envi to tif
            output format: *.tif
        3) Crop flightlines to AOP tile bounds
            output format: *.tif
        4) Merge tiffs that overlap spatially
            output format: *.tif
    """
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

    flightline_h5_path = os.path.join(data_raw_aop_path,site,year_aop,"hs_flightline")
    
    # 1) Apply BRDF/topographic corrections
    log.info(f'Applying BRDF/topo corrections for: {site} {year_inv}')
    # move_red_yellow_subset(site, flightline_h5_path) 
    # flightline_envi_path = convert_hdf_to_envi(flightline_h5_path, 
    #                                            os.path.join(data_int_path,site,year_inv,
    #                                                         'hs_flightline/'))
    flightline_envi_path = os.path.join(data_int_path,site,year_inv, 'hs_flightline/')
    # config_path = create_config(flightline_envi_path)
    # implement_brdf_correction(config_path) 

    # 2)  Convert corrected envi to tif
    log.info(f'Converting from envi to tif format for: {site} {year_inv}')
    convert_envi_to_tif(site, flightline_h5_path, flightline_envi_path) 

    # 3) Crop flightlines to AOP tile bounds
    log.info(f'Cropping flightline tifs to AOP tile bounds for: {site} {year_inv}')
    crop_flightlines_to_tiles = ro.r('crop_flightlines_to_tiles')
    # ???? = crop_flightlines_to_tiles(site, year_inv, data_raw_aop_path, flightline_envi_path) #ais

    # 4) Merge tiffs that overlap spatially
    log.info(f'Merging spatially overlapping tiffs for: {site} {year_inv}')
    

def prep_manual_training_data(site, year, data_raw_inv_path, data_int_path, biomass_path):
    """
    Clean and organize manually delineated tree crowns with PFT labels
    """

    # Create features (points or polygons) for each tree 
    log.info(f'Creating tree crown training data features for: {site} {year}')
    r_source = ro.r['source']
    # r_source(str(Path(__file__).resolve().parent/'inventory_helper.R'))
    # match_species_to_pft = ro.r('match_species_to_pft') 
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))
    # list_tiles_w_veg = ro.r('list_tiles_w_veg') 
    
    # Create tree polygons
    prep_manual_crown_delineations = ro.r('prep_manual_crown_delineations') 
    training_shp_path = prep_manual_crown_delineations(site, year, data_raw_inv_path, data_int_path, 
                                                   biomass_path)
    
    #ais abandoned converting this function from R to python because loading shapefiles with gpd lost the column data types. also tough to use list_tiles_w_veg function in python
    
    # training_data_dir = os.path.join(data_int_path, site, year, "training")
    # if not os.path.exists(training_data_dir):
    #     os.makedirs(training_data_dir)

    # # Merge manual shapes with NEON metadata for each site-year 
    # # Here we get the right PFT for each shape in each year

    # # NEON inventory for SJER 2017, SJER 2021, SOAP 2019, SOAP 2021, TEAK 2021                                    
    # neon_inv_manual = gpd.read_file(os.path.join(data_raw_inv_path,"manually_uploaded", f"NEON - {site} - {year}", "my_shapes.shp"))

    # # Load cleaned inventory data for site and year
    # # ais this is live trees only. should figure out how to incorporate dead trees
    # veg_ind = pd.read_csv(os.path.join(biomass_path, "pp_veg_structure_IND_IBA_IAGB_live.csv")).rename(columns={
    #     'individualID': 'indvdID',
    #     'adjEasting': 'easting',
    #     'adjNorthing': 'northing',
    #     'adjDecimalLatitude': 'lat',
    #     'adjDecimalLongitude': 'long',
    #     'adjElevation': 'elev',
    #     'individualStemNumberDensity': 'indStemDens',
    #     'individualBasalArea': 'indBA',
    #     'basalStemDiameter': 'bStemDiam',
    #     'basalStemDiameterMsrmntHeight': 'bStemMeasHgt',
    #     'scientificName': 'sciNameFull',
    #     'scientific': 'sp_gen',
    #     'wood.dens': 'wood_dens1',
    #     'wood.dens': 'wood_dens2'
    # }).assign(pft=lambda df: df['taxonID'].apply(match_species_to_pft))

    # # link inventory csv and manual shapes
    # neon_inv = neon_inv_manual.merge(veg_ind, on='indvdID', how='left')
    # neon_inv['id'] = neon_inv['id'].astype(str)

    # if site == "SJER" and year == "2021": 
    #     SJER_2021_neon_ref = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "NEON - SJER - 2021", "tree_crowns_training.shp"))
        
    #     # SJER 2023 carissa/anna fieldwork
    #     SJER_2023_c_manual = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "Carissa - SJER - 2023", "my_shapes.shp"))
    #     SJER_2023_c_ref = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "Carissa - SJER - 2023", "tgis_merged_species_32611.shp"))
        
    #     # SJER 2024 carissa/anna fieldwork
    #     SJER_2024_c_a_manual = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "Carissa Anna - SJER - 2024", "my_shapes.shp"))
    #     SJER_2024_c_a_ref1 = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "Carissa Anna - SJER - 2024", "trimble_merged_PLY_SJER.shp")).to_crs(epsg=32611)
    #     # SJER_2024_a_ref2 = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "Carissa Anna - SJER - 2024", "Anna_V1_ais-point.shp"))
    #     # SJER_2024_c_a_ref3 = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "Carissa Anna - SJER - 2024", "garmin_merged_PT_SJER.shp")).to_crs(epsg=32611)

    #     SJER_2021_neon = neon_inv.drop(columns=['id']).merge(SJER_2021_neon_ref[['indvdID', 'sp_gen']], on='indvdID').rename(columns={'sp_gen': 'species'})
    #     SJER_2023_c = SJER_2023_c_manual.merge(SJER_2023_c_ref[['Name', 'species']], on='Name').rename(columns={'Name': 'indvdID'})
    #     SJER_2024_c_a = SJER_2024_c_a_manual.merge(SJER_2024_c_a_ref1[['Name', 'Species', 'If other s', 'Site type']], on='Name').assign(
    #         species=lambda df: df.apply(lambda row: row['species_ai'] if pd.notna(row['species_ai']) else (
    #             'grass' if row['Site type'] == 'Grassland' else (
    #                 row['If other s'] if row['Species'] == 'OTHER' else row['Species'])), axis=1)).fillna({'species': 'grass'}).reset_index().assign(
    #         Name=lambda df: df.apply(lambda row: f"{row.name}_SJER24_ais" if pd.isna(row['Name']) else row['Name'], axis=1)).rename(columns={'Name': 'indvdID'})
        
    #     # Merge shapefiles for all years since we're using only 2021 imagery for now
    #     my_shapes = pd.concat([SJER_2021_neon, SJER_2023_c, SJER_2024_c_a]).assign(pft=lambda df: df.apply(lambda row: match_species_to_pft(row['species'], row['indvdID']), axis=1)).query("indvdID != 'A08'")

    # elif site == "TEAK" and year == "2021":
    #     TEAK_2021_neon_ref = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "NEON - TEAK - 2021", "tree_crowns_training.shp"))
        
    #     # TEAK 2023 carissa/anna fieldwork
    #     TEAK_2023_c_a_manual = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "Carissa - TEAK - 2023", "my_shapes.shp"))
    #     TEAK_2023_c_a_ref = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "Carissa - TEAK - 2023", "tgis_merged_species.shp")).to_crs(epsg=32611)
        
    #     # TEAK 2024 carissa/anna fieldwork
    #     TEAK_2024_c_a_manual = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "Carissa - TEAK - 2024", "my_shapes.shp"))
    #     TEAK_2024_c_a_ref = gpd.read_file(os.path.join(data_raw_inv_path, "manually_uploaded", "Carissa - TEAK - 2024", "tgis_merged_PLY_TEAK.shp")).to_crs(epsg=32611)

    #     TEAK_2021_neon = neon_inv.drop(columns=['id', 'notes']).merge(TEAK_2021_neon_ref[['indvdID', 'sp_gen']], on='indvdID').rename(columns={'sp_gen': 'species'})
        
    #     TEAK_2023_c_a = TEAK_2023_c_a_manual.merge(TEAK_2023_c_a_ref[['Name', 'species']], left_on='indvdID', right_on='Name').assign(
    #         species=lambda df: df.apply(lambda row: 'oak' if 'C03' in row['Name'] else (
    #             'willow' if 'C036' in row['Name'] else (
    #                 'pine' if 'A172' in row['Name'] else (
    #                     'manzanita' if 'A101' in row['Name'] else row['species']))), axis=1)).query("species.notna()").query("species.str.contains('lotus', case=False) == False").rename(columns={'Name': 'indvdID'})
        
    #     TEAK_2024_c_a = TEAK_2024_c_a_manual.rename(columns={'species_ai': 'species'}).reset_index().assign(
    #         indvdID=lambda df: df.apply(lambda row: f"{row.name}_TEAK24_ais" if pd.isna(row['id']) else row['id'], axis=1)).query("species.notna()")

    #     # Merge shapefiles for all years since we're using only 2021 imagery for now
    #     my_shapes = pd.concat([TEAK_2021_neon, TEAK_2023_c_a, TEAK_2024_c_a]).assign(pft=lambda df: df.apply(lambda row: match_species_to_pft(row['species'], row['indvdID']), axis=1)).drop(columns=['geometry'])

    # else:
    #     # rename df so that list_tiles_w_veg is generalizable below
    #     my_shapes = neon_inv

    # # Generate a list of AOP tiles that overlap with the inventory data
    # # my_shapes_r = pandas2ri.py2rpy(my_shapes)
    # # list_tiles_w_veg(veg_df = my_shapes_r, out_dir = training_data_dir)


    # my_shapes['domainID'] = my_shapes['domainID'].astype(str)
    # my_shapes['siteID'] = my_shapes['siteID'].astype(str)
    # my_shapes['species_ai'] = my_shapes['species_ai'].astype(str)
    # my_shapes['Species'] = my_shapes['Species'].astype(str)
    # my_shapes['If other s'] = my_shapes['If other s'].astype(str)

    # training_shp_path = os.path.join(training_data_dir, "ref_labelled_crowns.shp")
    # my_shapes.to_file(training_shp_path)

    log.info('Clipped tree crown polygons data saved at: '
             f'{training_shp_path}')
    
    return training_shp_path


def prep_aop_imagery(site, year, hs_type, hs_path, tif_path, data_int_path, use_tiles_w_veg):
    """ Prepare aop imagery for extracting descriptive features for classification by 
            (1) clean/mask out shadow and non-veg
            (2) stacking imagery 
     Analog to 05-prep_aop_imagery.R from https://github.com/earthlab/neon-veg 

     ais to do: Then plots one tile of the prepped imagery. Analog to 06-plot_aop_imagery.R

    Calls a function written in R

    Parameters
    ----------
    site : str
        Site name
    year : str
        Year of aop collection
    hs_path : str
        Path where either raw hyperspectral h5 files OR BRDF-corrected tifs are stored
    tif_path : str
        Path where the raw aop tif files are stored

    Returns
    -------
    (str)
        Path to the result folder for stacked HS data ready for classification
    """
    log.info(f'Prepping aop data for: {site} {year}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

    # Clean/ mask imagery

    # Stack imagery
    prep_aop_imagery = ro.r('prep_aop_imagery')
    stacked_aop_path = prep_aop_imagery(site, year, hs_type, hs_path, tif_path, data_int_path, 
                                        use_tiles_w_veg)
    log.info('Stacked AOP data and saved at: '
             f'{stacked_aop_path}')
    return stacked_aop_path


def extract_spectra_from_polygon(site, year, shp_path, data_int_path, data_final_path, stacked_aop_path, 
                         use_case, aggregate_from_1m_to_2m_res, ic_type):
    """Create geospatial features (points, polygons with half the maximum crown diameter) 
        for every tree in the NEON woody vegetation data set. 
        Analog to 02-create_tree_features.R from https://github.com/earthlab/neon-veg 

        Generate polygons that intersect with independent pixels in the AOP data 
        Analog to 03-process_tree_features.R from https://github.com/earthlab/neon-veg 

        Extract features (remote sensing data) for each sample (pixel) within the 
        specified shapefile (containing polygons that correspond to trees at the NEON site)
        Analog to 07-extract_training_features.R from https://github.com/earthlab/neon-veg 
    """

    # # Create features (points or polygons) for each tree 
    # log.info(f'Extracting spectra from tree crowns for: {site} {year}')

    # training_data_dir = os.path.join(data_int_path, site, year, "training")
    # ic_type_path = os.path.join(data_final_path,site,year,ic_type)

    # # get a description of the shapefile to use for naming outputs
    # shapefile_description = os.path.splitext(os.path.basename(shp_path))[0]

    # # Specify destination for extracted features
    # if use_case == "train":
    #     extracted_features_path     = os.path.join(training_data_dir)
    #     extracted_features_filename = os.path.join(extracted_features_path,
    #                                              shapefile_description+"-extracted_features_inv.tif")                
    # elif use_case=="predict":
    #     extracted_features_path     = os.path.join(ic_type_path) 
    #     extracted_features_filename = os.path.join(extracted_features_path,
    #                                              shapefile_description+"-extracted_features.tif")
    # else:
    #     print("need to specify use_case")

    # # Only run if the extracted features do not exist
    # if not os.path.exists(extracted_features_filename):
    #     # Load shapefile
    #     shp_gdf = gpd.read_file(shp_path)

    #     # Compute centroid for each geometry
    #     shp_gdf["center_X"] = shp_gdf.geometry.centroid.x
    #     shp_gdf["center_Y"] = shp_gdf.geometry.centroid.y

    #     # List all .tif raster files in the directory
    #     stacked_aop_list = [os.path.join(stacked_aop_path, f) for f in os.listdir(stacked_aop_path) if f.endswith(".tif")]

    #     # create column to track shape ID 
    #     # if training, this is the tree crown boundary
    #     # if predicting, this is the plot boundary)
    #     if use_case == "train":
    #         shp_gdf["shapeID"] = [f"tree_crown_{i}" for i in range(len(shp_gdf))]
    #     elif use_case == "predict":
    #         if "PLOTID" in shp_gdf.columns:
    #             shp_gdf = shp_gdf.rename(columns={"PLOTID": "shapeID"})
    #         elif "plotID" in shp_gdf.columns:
    #             shp_gdf = shp_gdf.rename(columns={"plotID": "shapeID"})
    #         else:
    #             print("Trying to predict without training data. Need to first switch use_case to train from predict")

    #     # Filter to tiles containing veg to speed up the next for-loop
    #     # ais does this code include tiles with only a crown from the ucla field data though?
    #     if use_case == 'train' or ic_type=='rs_inv_plots':
    #         # Read tiles_w_veg.txt and extract the first column as a list
    #         tiles_w_veg = pd.read_csv(os.path.join(training_data_dir, "tiles_w_veg.txt"), header=None)[0].tolist()
                    
    #         # Filter stacked_aop_list based on the presence of any tile name
    #         stacked_aop_list = [path for path in stacked_aop_list if any(tile in path for tile in tiles_w_veg)]

    #     # Loop through AOP tiles
    #     for stacked_aop_filename in stacked_aop_list:

    #         # Read current tile of stacked AOP data
    #         with rasterio.open(stacked_aop_filename) as src:
    #             stacked_aop_data = src.read()
    #             raster_transform = src.transform
    #             raster_crs = src.crs
    #             raster_bounds = src.bounds  # Get raster bounds

    #             # Construct the easting northing string for naming outputs
    #             east_north_string = f"{round(raster_bounds[0])}_{round(raster_bounds[1])}"
                
    #             east_north_tif_path = os.path.join(extracted_features_path, f"extracted_features_mask_{east_north_string}_{shapefile_description}.tif")
                
    #             # If CSV file already exists, skip
    #             if os.path.exists(east_north_tif_path):
    #                 continue

    #             # Check which polygons overlap with the raster
    #             shp_gdf = shp_gdf[shp_gdf.intersects(box(*raster_bounds))]

    #             # If no polygons overlap, exit
    #             if shp_gdf.empty:
    #                 print("No overlapping polygons found.")
    #                 continue
    #             else:
    #                 # Convert overlapping polygons to GeoJSON format
    #                 shapes_geojson = [mapping(geom) for geom in shp_gdf.geometry]
                    
    #                 # Mask the raster using the overlapping polygons
    #                 with rasterio.open(stacked_aop_filename) as src:
    #                     masked_raster, masked_transform = rasterio.mask.mask(dataset=src, shapes=shapes_geojson, crop=False, nodata=np.nan)
                    
    #                 # Save the masked raster as a new GeoTIFF file
    #                 with rasterio.open(
    #                     east_north_tif_path, "w", driver="GTiff",
    #                     height=masked_raster.shape[1], width=masked_raster.shape[2],
    #                     count=masked_raster.shape[0], dtype=str(masked_raster.dtype),
    #                     crs=raster_crs, transform=masked_transform
    #                 ) as dst:
    #                     dst.write(masked_raster)
                    
    #                 print(f"Masked raster saved to {east_north_tif_path}")

    #     # combine all extracted features into a single .tif
    #     paths_ls = glob.glob(os.path.join(extracted_features_path, "*.tif"))
        
    #     # refine the output csv selection 
    #     tifs = [path for path in paths_ls if f"000_{shapefile_description}.tif" in path]
        
    #     # Open and merge TIF files
    #     src_files_to_mosaic = [rasterio.open(tif) for tif in tifs]
    #     mosaic, out_trans = merge(src_files_to_mosaic)
    #     # Open and merge TIF files
    #     src_files_to_mosaic = [rasterio.open(tif) for tif in tifs]
    #     mosaic, out_trans = merge(src_files_to_mosaic)

    #     # Get metadata from the first file
    #     out_meta = src_files_to_mosaic[0].meta.copy()   
    #     # Get metadata from the first file
    #     out_meta = src_files_to_mosaic[0].meta.copy()   

    #     # Update metadata for the merged file
    #     out_meta.update({
    #         "driver": "GTiff",
    #         "height": mosaic.shape[1],
    #         "width": mosaic.shape[2],
    #         "transform": out_trans
    #     })

    #     # Write the merged TIF file
    #     with rasterio.open(extracted_features_filename, "w", **out_meta) as dest:
    #         dest.write(mosaic)

    #     # Close input files
    #     for src in src_files_to_mosaic:
    #         src.close()

    #     # Delete the individual TIF files for each tile
    #     for tif in tifs:
    #         os.remove(tif)

    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))
    extract_spectra_from_polygon_r = ro.r('extract_spectra_from_polygon_r')

    # Extract training data from AOP data with tree polygons
    training_spectra_path = extract_spectra_from_polygon_r(site=site, 
                                                         year=year, 
                                                         data_int_path=data_int_path, 
                                                         data_final_path=data_final_path, 
                                                         stacked_aop_path=stacked_aop_path, 
                                                         shp_path=shp_path,
                                                         use_case=use_case, 
                                                         aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res,
                                                         ic_type=ic_type)

    log.info('Spectral features for training data saved at: '
             f'{training_spectra_path}')
    
    return training_spectra_path



def extract_spectra_3darray(raster_path, training_shp):
    # generated by cborg, adapted by AIS March 2025

    # Load the multiband raster
    with rasterio.open(raster_path) as src:
        # Prepare a 3D array to store training data (bands, height, width)
        training_array = []

        # Create a combined mask for all polygons
        combined_mask = np.zeros(src.shape, dtype=bool)

        for _, row in training_shp.iterrows():
            geometry = row['geometry']
            # Create a mask for the raster using the polygon
            mask = geometry_mask([geometry], invert=True, transform=src.transform, out_shape=src.shape)

            # Combine the mask
            combined_mask |= mask  # Logical OR to combine masks

        # Read the raster values
        raster_values = src.read()

        # Extract pixel values for the masked area
        for band in range(raster_values.shape[0]):
            band_values = raster_values[band][combined_mask]
            training_array.append(band_values)

        # Convert the list of arrays into a 3D numpy array
        training_array = np.array(training_array)  # Shape will be (bands, num_samples)

        return training_array, combined_mask



def train_pft_classifier(sites, data_int_path, pcaInsteadOfWavelengths, ntree, 
                         randomMinSamples, independentValidationSet):    
    """Train a Random Forest (RF) model to classify tree PFT using in-situ tree
        measurements for PFT labels and remote sensing data as descriptive features
        Analog to 08-classify_species.R from https://github.com/earthlab/neon-veg 

        Train random forest model and assess PFT classification accuracy.
        Model outputs will be written to a folder within the training directory 
        starting with "rf_" followed by a description of each shapefile 
        containing points or polygons per tree. 

        Evaluate classifier performance: Out-of-Bag accuracy, independent
        validation set accuracy, and Cohen's kappa. Generate confusion matrix
        to show the classification accuracy for each PFT. 
        Analog to 09-assess_accuracy.R from https://github.com/earthlab/neon-veg 

        ntree 
          RF tuning parameter, number of trees to grow. default value 500
    
        randomMinSamples 
          To reduce bias, this boolean variable indicates whether the same number 
          of samples are randomly selected per PFT class. Otherwise,
          all samples per class are used for training the classifier. 
    
        independentValidationSet=T
          if TRUE, keep separate set for validation
          ais this doesnt work if for F - troubleshoot this later

    """
    
    log.info(f'Training model on sites: {sites}')

    # Directory to save the classification outputs 
    rf_output_dir = os.path.join(data_int_path, "rf_dir")
    if not os.path.exists(rf_output_dir):
        os.makedirs(rf_output_dir) 

    features_df = pd.DataFrame()
    for site in sites:
        site_year_paths = [os.path.join(data_int_path, site, year) for year in os.listdir(os.path.join(data_int_path, site))]
        for site_year_path in site_year_paths:

            if not os.path.exists(os.path.join(site_year_path, "stacked_aop/wavelengths.txt")):
                continue

            # Load data
            wavelengths = pd.read_csv(os.path.join(site_year_path, "stacked_aop/wavelengths.txt"))
            # Stacked AOP layer names
            stacked_aop_layer_names = pd.read_csv(os.path.join(site_year_path, "stacked_aop/stacked_aop_layer_names.txt"))
            # Labelled, half-diam crown polygons to be used for training/validation
            shapefile_description = os.path.basename(os.path.join(site_year_path, "training/ref_labelled_crowns.shp")).split('.')[0]
            # Csv file containing extracted features
            extracted_features_filename = os.path.join(site_year_path, "training/ref_labelled_crowns-extracted_features_inv.csv")
            if not os.path.exists(extracted_features_filename):
                continue

            # Filter out unwanted wavelengths
            wavelength_lut = filter_out_wavelengths(wavelengths=wavelengths['wavelengths'].tolist(), 
                                                    layer_names=stacked_aop_layer_names['stacked_aop_layer_names'].tolist())
            
            # features and label to use in the RF models
            featureNames = ["shapeID", "pft"] + wavelength_lut['xwavelength'].tolist() + [name for name in stacked_aop_layer_names['stacked_aop_layer_names'].tolist() if not name.isdigit()]
        
            # Prep extracted features csv for RF model
            extracted_features_X_df = pd.read_csv(extracted_features_filename)
            extracted_features_X_df = extracted_features_X_df[featureNames] 
            # Remove X from wavelength column names
            features_temp = extracted_features_X_df.rename(columns=lambda x: x[1:] if x.startswith('X') else x)                  
            # filter the data to contain only the features of interest 
            # features_temp = extracted_features_df[featureNames] 
            # features_temp.columns = featureNames[len(wavelengths):] + [f"X{i+1}" for i in range(len(wavelengths))] 
            
            # make sure that column names are the same so we can rbind them across site-year
            numeric_columns = [col for col in features_temp.columns if col.isdigit()]
            features_temp.columns = features_temp.columns.map(lambda x: f'X{numeric_columns.index(x)+1}' if x in numeric_columns else x)

            features_df = pd.concat([features_df, features_temp])

    # Remove any rows with NA   
    features_df.dropna(inplace=True)
        
    # group the data by 'pft' and 'shapeID'
    grouped_features = features_df.groupby(['pft', 'shapeID'])
    groups = [group for _, group in grouped_features] # create a list of groups

    # split the groups into training and testing sets
    train_groups, test_groups = train_test_split(groups, test_size=0.2, random_state=42)
    train_df = pd.concat(train_groups)
    test_df = pd.concat(test_groups)

    # create the training and testing DataFrames
    train_df = train_df.drop(["shapeID"],axis=1) #'site', 'year'
    test_df  = test_df.drop(["shapeID"],axis=1) #'site','year'
    # check the class distribution in the training and testing sets
    print("Training set class distribution:")
    print(train_df['pft'].value_counts(normalize=True))
    print("Testing set class distribution:")
    print(test_df['pft'].value_counts(normalize=True))
    
    # Standardize the Data
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(train_df.loc[:, train_df.columns.str.startswith('X')])

    # Apply PCA
    pca = PCA(n_components=0.99)
    X_train_pca = pca.fit_transform(X_train_scaled) 

    # Save scaler and pca models to apply to new data for prediction
    joblib.dump(scaler, os.path.join(rf_output_dir,'scaler.joblib'))
    joblib.dump(pca, os.path.join(rf_output_dir,'pca.joblib'))

    # remove the individual spectral reflectance bands from the training data
    train_df_noWavelengths = train_df.drop([col for col in train_df.columns if col.startswith("X")], axis=1)#ais one time I won't need this line after I stop using 'stacekd_aop_layers'
    train_features = pd.concat([train_df_noWavelengths.reset_index(drop=True), pd.DataFrame(X_train_pca, columns=[f"PC{i+1}" for i in range(X_train_pca.shape[1])])],axis=1)
    
    if X_train_pca.shape[1] > 2:    
        # visualize PCA    
        plt.figure(figsize=(8, 6))
        for value in train_features['pft'].unique():
            plt.scatter(train_features.loc[train_features['pft'] == value, 'PC1'], train_features.loc[train_features['pft'] == value, 'PC2'], label=value)
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title('PC1 vs PC2')
        plt.legend()
        plt.savefig(os.path.join(rf_output_dir, 'pc1_vs_pc2.png')) # plot PC1 vs PC2
        
        plt.figure(figsize=(8, 6))
        for value in train_features['pft'].unique():
            plt.scatter(train_features.loc[train_features['pft'] == value, 'PC2'], train_features.loc[train_features['pft'] == value, 'PC3'], label=value)
        plt.xlabel('PC2')
        plt.ylabel('PC3')
        plt.title('PC2 vs PC3')
        plt.legend()
        plt.savefig(os.path.join(rf_output_dir, 'pc2_vs_pc3.png')) # plot PC2 vs PC3

        if X_train_pca.shape[1] > 3:       
            plt.figure(figsize=(8, 6))
            for value in train_features['pft'].unique():
                plt.scatter(train_features.loc[train_features['pft'] == value, 'PC3'], train_features.loc[train_features['pft'] == value, 'PC4'], label=value)
            plt.xlabel('PC3')
            plt.ylabel('PC4')
            plt.title('PC3 vs PC4')
            plt.legend()
            plt.savefig(os.path.join(rf_output_dir, 'pc3_vs_pc4.png')) 

    # create a bar plot to visualize the count of data points in training and testing datasets for each prediction class
    plt.figure(figsize=(10, 6))
    train_counts = train_df['pft'].value_counts()
    test_counts = test_df['pft'].value_counts()
    plt.bar(train_counts.index, train_counts.values, label='Training')
    plt.bar(test_counts.index, test_counts.values, label='Testing')
    plt.xlabel('Prediction Class')
    plt.ylabel('Count')
    plt.title('Count of Data Points in Training and Testing Datasets')
    plt.legend()
    plt.savefig(os.path.join(rf_output_dir, 'train_test_count.png'))

    if randomMinSamples:
        minSamples = train_df["pft"].value_counts().min()
        train_df = train_df.groupby("pft").apply(lambda x: x.sample(minSamples)).reset_index(drop=True)

    rf_model_path = os.path.join(rf_output_dir, f"rf_model.joblib")
    if not os.path.exists(rf_model_path):
        rf_model = RF(n_estimators=ntree, random_state=42, class_weight="balanced")
        rf_model.fit(train_features.drop(["pft"], axis=1), train_features["pft"])
        dump(rf_model, rf_model_path)
        log.info('Trained PFT classifier saved in this folder: '
                f'{rf_model_path}')
    else:
        rf_model = load(rf_model_path)

    # Prepare test data
    X_test_scaled = scaler.transform(test_df.loc[:, test_df.columns.str.startswith('X')])  # Use the same scaler
    X_test_pca    = pca.transform(X_test_scaled)    
    # remove the individual spectral reflectance bands from the training data
    test_df_noWavelengths = test_df.drop([col for col in test_df.columns if col.startswith("X")], axis=1)#ais one time I won't need this line after I stop using 'stacekd_aop_layers'
    test_features = pd.concat([test_df_noWavelengths.reset_index(drop=True), pd.DataFrame(X_test_pca, columns=[f"PC{i+1}" for i in range(X_test_pca.shape[1])])],axis=1)
    
    # Predict test data
    y_pred = rf_model.predict(test_features.drop(["pft"], axis=1))

    # Calculate metrics
    conf_matrix = confusion_matrix(test_features["pft"], y_pred)
    cm_df = pd.DataFrame(conf_matrix, index=sorted(set(y_pred)), columns=sorted(set(y_pred)))
    accuracy = accuracy_score(test_features["pft"], y_pred)
    precision = precision_score(test_features["pft"], y_pred, average='weighted')
    recall = recall_score(test_features["pft"], y_pred, average='weighted')
    f1 = f1_score(test_features["pft"], y_pred, average='weighted')

    summary_stats = [
        "Confusion Matrix:\n" + cm_df.to_string(),  # Convert DataFrame to string
        f"\nAccuracy: {accuracy:.4f}",
        f"Precision: {precision:.4f}",
        f"Recall: {recall:.4f}",
        f"F1 Score: {f1:.4f}",
        "\nClass-level Report:\n" + classification_report(test_features["pft"], y_pred),
        ]

    # Write summary statistics to the file
    rf_summ_txt_path = os.path.join(rf_output_dir,'rf_summary_statistics.txt')
    with open(rf_summ_txt_path, "w") as f:
            f.write("\n".join(summary_stats))  
        
    # Print to terminal window
    print("Validation Accuracy:", accuracy)
    print("Classification Report:")
    print(classification_report(test_features["pft"], y_pred))
    print("Confusion Matrix:")
    print(cm_df)

    # Plot variable importance
    importances = rf_model.feature_importances_
    importance_df_unsorted = pd.DataFrame({"Feature": test_features.drop(["pft"], axis=1).columns, "Importance": importances})
    importance_df = importance_df_unsorted.sort_values(by="Importance", ascending=False)

    # Plot Variable Importance
    var_imp_path = os.path.join(rf_output_dir, "feature_importance_plot.png")
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Importance", y="Feature", data=importance_df, palette="viridis")
    plt.xlabel("Feature Importance Score")
    plt.ylabel("Features")
    plt.title("Feature Importance from Random Forest")
    plt.savefig(var_imp_path, dpi=300, bbox_inches="tight")  # High-resolution save
    
    return rf_model_path   



def plot_cv_indices(cv, X, y, groups, ax, n_splits, cmap_data, cmap_cv, lw=10):
    """Create a sample plot for indices of a cross-validation object.
    from: https://scikit-learn.org/stable/auto_examples/model_selection/plot_cv_indices.html"""
    #use_groups = "Group" in type(cv).__name__
    # groups = group if use_groups else None
    # Generate the training/testing visualizations for each CV split
    for ii, (tr, tt) in enumerate(cv.split(X=X, y=y, groups=groups)):
        # Fill in indices with the training/test groups
        indices = np.array([np.nan] * len(X))
        indices[tt] = 1
        indices[tr] = 0

        # Visualize the results
        ax.scatter(
            range(len(indices)),
            [ii + 0.5] * len(indices), 
            c=indices,
            marker="_",
            lw=lw,
            cmap=cmap_cv,
            vmin=-0.2,
            vmax=1.2,
        )

    # Plot the data classes and groups at the end
    ax.scatter(
        range(len(X)), [ii + 1.5] * len(X), c=y, marker="_", lw=lw, cmap=cmap_data
    )
    ax.scatter(
        range(len(X)), [ii + 2.5] * len(X), c=groups, marker="_", lw=lw, cmap=cmap_data
    )

    # Formatting
    yticklabels = list(range(n_splits)) + ["class", "group"]
    ax.set(
        yticks=np.arange(n_splits + 2) + 0.5,
        yticklabels=yticklabels,
        xlabel="Sample index",
        ylabel="CV iteration",
        ylim=[n_splits + 2.2, -0.2],
        # xlim=[0, 100],
    )
    ax.set_title("{}".format(type(cv).__name__), fontsize=15)
    return ax



def fit_RF_CV_class(X, y, groups, k_fold, pca, savefile):
    # Written by Nicola Falco
    # Adapted by Anna Spiers Nov 2024
  
    ##################### Define evaluation procedure for CV
    cv = StratifiedKFold(n_splits=k_fold)#, shuffle=True, random_state=10210)

    ##################### test/training 
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=111)

    #############################################################################
    #####################  cross-validation
    if pca == True:
        pipeline = Pipeline(steps=[
            ['pca', PCA(n_components=0.99, whiten= True)],
            # ['scaler', MinMaxScaler()],
            ['classifier', RF(class_weight='balanced')] 
            ])
        savefile = (str(savefile) + '_PCA')

    else:
        pipeline = Pipeline(steps=[
            # ['scaler', MinMaxScaler()],
            ['classifier', RF(class_weight='balanced')]
            ])
        savefile = (str(savefile))
        
    ##################### Define grid search RF Parameter optimization
    # Number of trees in random forest
    n_estimators = [int(x) for x in np.linspace(start = 1000, stop = 10000, num = 4)] #5

    # Number of features to consider at every split
    max_features = ['auto', 'sqrt']

    # Maximum number of levels in tree
    max_depth = [int(x) for x in np.linspace(10, 110, num = 4)] #5
    max_depth.append(None)

    # Minimum number of samples required to split a node
    # min_samples_split = [2, 5, 10]

    # Minimum number of samples required at each leaf node
    min_samples_leaf = [1, 2, 4]

    # Method of selecting samples for training each tree
    bootstrap = [True]#, False]

    param_rf ={'classifier__n_estimators': n_estimators,
               'classifier__max_features': max_features,
               'classifier__min_samples_leaf': min_samples_leaf,
               'classifier__max_depth': max_depth,
               # 'classifier__min_samples_split': min_samples_split,
               # 'classifier__min_samples_leaf': min_samples_leaf,
               'classifier__bootstrap': bootstrap
               }

    search = GridSearchCV(estimator=pipeline, param_grid=param_rf,
                           scoring='f1_weighted', cv=cv, refit=True, verbose=3, n_jobs = 128)
  
    ##################### Fit the search to the data
    # Fit the grid search to the data
    search.fit(X_train, y_train)
    search.best_params_

    ##################### get the best performing model fit on the whole training set
    best_model = search.best_estimator_
  
    ##################### Prediction
    p_test = best_model.predict(X_test)

    ##################### Report performance
    #### average CM and classification report
    class_list= np.unique(y_test)

    cm_analysis(y_test, p_test, class_list, savefile)
    class_report(y_test, p_test, class_list,savefile)
    
    return best_model, X_train, X_test



def cm_analysis(y_true, y_pred, labels,savefile): 

    SMALL_SIZE = 12
    MEDIUM_SIZE = 15
    BIGGER_SIZE = 18
   
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
 
    sns.set(style="whitegrid")

    ##############################
 
    # compute the CM
    cm = confusion_matrix(y_true, y_pred, labels=labels, normalize='true')
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=labels)
    fig, ax = plt.subplots(figsize=(8,8))
    disp.plot(ax=ax, cmap=plt.cm.Blues)

    # save the plot as PNG
    save2png= (str(savefile) + '_CMnorm.png')
    plt.savefig(save2png, dpi=300, bbox_inches='tight')    

    # transform cm into dataframe and save and CSV
    df_cm = pd.DataFrame(cm)
    save2csv= (str(savefile) + '_CMnorm.csv')
    df_cm.to_csv(save2csv)

    

def class_report(y_test, y_pred, labels,savefile):
    report = classification_report(y_test,y_pred,labels=labels,output_dict=True)
    df_report = pd.DataFrame(report).transpose()
    save2csv = (str(savefile) + '_ClassReport.csv')
    df_report.to_csv(save2csv)




def serial_index_to_coordinates(serial_index, raster_width=1000):
    """
    Translates a serial index of a raster into coordinates of the raster.

    Parameters:
    serial_index (int): The serial index of the pixel in the raster.
    raster_width (int): The width of the raster.

    Returns:
    tuple: A tuple containing the x and y coordinates of the pixel in the raster.
    """
    x = (serial_index // raster_width)-1
    y = 1000 - (serial_index % raster_width) 
    return x, y





def filter_out_wavelengths(wavelengths, layer_names):
    # define the "bad bands" wavelength ranges in nanometers, where atmospheric 
    # absorption creates unreliable reflectance values. 
    # bad_band_window_1 = (1340, 1445)
    # bad_band_window_2 = (1790, 1955)
    wavelengths_np = np.array(wavelengths)

    # remove the bad bands from the list of wavelengths 
    remove_bands = wavelengths_np[(wavelengths_np > 1340) & 
                               (wavelengths_np < 1445) | 
                               (wavelengths_np > 1790) & 
                               (wavelengths_np < 1955)]
    
    # Make sure printed wavelengths and stacked AOP wavelengths match
    if not np.allclose(np.round(wavelengths_np), np.array([float(name) for name in layer_names[:len(wavelengths_np)]])):
        print("wavelengths do not match between wavelength.txt and the stacked imagery")
    
    # create a LUT that matches actual wavelength values with the column names,
    # X followed by the rounded wavelength values. 
    # Remove the rows that are within the bad band ranges. 
    wavelength_lut = pd.DataFrame({'wavelength': wavelengths_np,
                                  'xwavelength': ['X' + str(round(w)) for w in wavelengths_np]})[~np.in1d(wavelengths_np, remove_bands)]
    
    return wavelength_lut


