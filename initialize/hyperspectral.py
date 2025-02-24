import logging
import os
import sys
# import whitebox
import zipfile
import ray

import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Patch
import rasterio
from rasterio.plot import show
import fiona
from sklearn.preprocessing import OneHotEncoder
# from keras.models import Sequential
# from keras.layers import Conv2D, MaxPooling2D, Flatten, Dense
import geopandas as gpd
from shapely.geometry import Point
import rpy2.robjects as ro
from collections import Counter
import seaborn as sns

# data preparation
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

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

    product_codes = ['DP3.30015.001', 'DP3.30010.001'] #chm, rgb
    file_type = 'tif' 
    for product_code in product_codes: 
        p = path/file_type
        download_aop_files(product_code,
                           site,
                           year,
                           str(p),
                           match_string=file_type,
                           check_size=False)

    product_code = 'DP3.30026.001' #veg indices
    p = path/file_type
    download_aop_files(product_code,
                       site,
                       year,
                       str(p),
                       match_string='.zip',
                       check_size=False)
    file_type = 'tif'
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

def generate_pft_reference(sites, data_raw_inv_path, data_int_path, trait_table_path):
    """ Generate a csv that translates each species into a PFT usable as training data
    """
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

    # Generate PFT reference
    generate_pft_reference = ro.r('generate_pft_reference')
    Path(data_int_path).mkdir(parents=True, exist_ok=True)
    pft_reference_path = generate_pft_reference(sites, data_raw_inv_path, data_int_path, trait_table_path)
    log.info('Generated PFT reference for: {sites}'
             f'{pft_reference_path}')
    #ais now where to feed in pft_reference csv path


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
    

def prep_manual_training_data(site, year, data_raw_inv_path, data_int_path, biomass_path, 
                               pft_reference_path, px_thresh):
    """
    Clean and organize manually delineated tree crowns with PFT labels
    """

    # Create features (points or polygons) for each tree 
    log.info(f'Creating tree crown training data features for: {site} {year}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))
    
    # Create tree polygons
    prep_manual_crown_delineations = ro.r('prep_manual_crown_delineations') 
    training_shp_path = prep_manual_crown_delineations(site, year, data_raw_inv_path, data_int_path, 
                                                   biomass_path, pft_reference_path, px_thresh)
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
                         use_case, aggregate_from_1m_to_2m_res, ic_type=None, ic_type_path=None):
    """Create geospatial features (points, polygons with half the maximum crown diameter) 
        for every tree in the NEON woody vegetation data set. 
        Analog to 02-create_tree_features.R from https://github.com/earthlab/neon-veg 

        Generate polygons that intersect with independent pixels in the AOP data 
        Analog to 03-process_tree_features.R from https://github.com/earthlab/neon-veg 

        Extract features (remote sensing data) for each sample (pixel) within the 
        specified shapefile (containing polygons that correspond to trees at the NEON site)
        Analog to 07-extract_training_features.R from https://github.com/earthlab/neon-veg 
    """

    # Create features (points or polygons) for each tree 
    log.info(f'Extracting spectra from tree crowns for: {site} {year}')
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
                                                         aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res)
 

    log.info('Spectral features for training data saved at: '
             f'{training_spectra_path}')
    
    return training_spectra_path



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
            shapefile_description = os.path.basename(os.path.join(site_year_path, "training/tree_crowns_training.shp")).split('.')[0]
            # Csv file containing extracted features
            extracted_features_filename = os.path.join(site_year_path, "training/tree_crowns_training-extracted_features_inv.csv")
            if not os.path.exists(extracted_features_filename):
                continue

            # Directory to save the classification outputs 
            rf_output_dir = os.path.join(data_int_path, "rf_dir")
            if not os.path.exists(rf_output_dir):
                os.makedirs(rf_output_dir)
                os.makedirs(os.path.join(rf_output_dir,"validation"))
            
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

    # remove the individual spectral reflectance bands from the training data
    features_noWavelengths = features_df.drop([col for col in features_df.columns if col.startswith("X")], axis=1)
    features_noWavelengths.drop(columns=['pixelNumber','eastingIDs','northingIDs'], inplace=True)
       
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features_df[[col for col in features_df.columns if col.startswith("X")]])
    pca = PCA(n_components=0.99) #figure out why so much variance is explained by PC1
    features_pca = pca.fit_transform(features_scaled) 
    nPCs = features_pca.shape[1] 
    features = pd.concat([features_noWavelengths.reset_index(drop=True), pd.DataFrame(features_pca, columns=[f"PC{i+1}" for i in range(nPCs)])],axis=1)
    if nPCs > 2:    
        # visualize PCA    
        plt.figure(figsize=(8, 6))
        for value in features['pft'].unique():
            plt.scatter(features.loc[features['pft'] == value, 'PC1'], features.loc[features['pft'] == value, 'PC2'], label=value)
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title('PC1 vs PC2')
        plt.legend()
        plt.savefig(os.path.join(rf_output_dir, 'pc1_vs_pc2.png')) # plot PC1 vs PC2
        
        plt.figure(figsize=(8, 6))
        for value in features['pft'].unique():
            plt.scatter(features.loc[features['pft'] == value, 'PC2'], features.loc[features['pft'] == value, 'PC3'], label=value)
        plt.xlabel('PC2')
        plt.ylabel('PC3')
        plt.title('PC2 vs PC3')
        plt.legend()
        plt.savefig(os.path.join(rf_output_dir, 'pc2_vs_pc3.png')) # plot PC2 vs PC3

        if nPCs > 3:       
            plt.figure(figsize=(8, 6))
            for value in features['pft'].unique():
                plt.scatter(features.loc[features['pft'] == value, 'PC3'], features.loc[features['pft'] == value, 'PC4'], label=value)
            plt.xlabel('PC3')
            plt.ylabel('PC4')
            plt.title('PC3 vs PC4')
            plt.legend()
            plt.savefig(os.path.join(rf_output_dir, 'pc3_vs_pc4.png')) 




    # # Define evaluation procedure for CV
    # rng = np.random.RandomState(1338)
    # cmap_data = plt.cm.Paired
    # cmap_cv = plt.cm.coolwarm
    # n_splits = 4 #ais make n_splits a global param 
    # # example: https://stackoverflow.com/questions/56872664/complex-dataset-split-stratifiedgroupshufflesplit

    # # Reorder features by pft
    # features_sort = features.sort_values(['pft', 'shapeID'])
    # X = features_sort.drop(columns=['pft','shapeID']) #list(range(len(features)))
    # y = [list(set(features_sort['pft'])).index(pft) for pft in features_sort['pft'] ] #use index number
    # groups = [list(set(features_sort['shapeID'])).index(shapeID) for shapeID in features_sort['shapeID']]  #use index number  
    # # for train, test in sgkf.split(X, y, groups=groups):
    # #     print("%s %s" % (train, test))
    # # ais need to print the index of each class - otherwise hard to interpret results of 0-10
    # best_model, X_train, X_test = fit_RF_CV_class(X=X, y=y, groups=groups, k_fold=n_splits, pca=False, 
    #             savefile=os.path.join(rf_output_dir,'pipeline_output'))
    
    # # split the groups into training and testing sets
    # train_df = pd.concat(X_train)
    # test_df = pd.concat(X_test)

    # # ais need to put this code ^ into a for-loop to choose best hyperparameters

    # fig, ax = plt.subplots(figsize=(6, 3))
    # plot_cv_indices(StratifiedGroupKFold(n_splits), X, y, groups, ax, n_splits, cmap_data, cmap_cv)
    # ax.legend(
    #     [Patch(color=cmap_cv(0.8)), Patch(color=cmap_cv(0.02))],
    #     ["Testing set", "Training set"],
    #     loc=(1.02, 0.8),
    # ) #ais adapt this code to add key for class
    # # Make the legend fit
    # plt.tight_layout()
    # fig.subplots_adjust(right=0.7)
    # plt.savefig(os.path.join(rf_output_dir,'StratGroupKFold_split.png'))




    # group the data by 'pft' and 'shapeID'
    grouped_features = features.groupby(['pft', 'shapeID'])
    groups = [group for _, group in grouped_features] # create a list of groups

    # split the groups into training and testing sets
    train_groups, test_groups = train_test_split(groups, test_size=0.2, random_state=42)
    train_df = pd.concat(train_groups)
    test_df = pd.concat(test_groups)




    # create the training and testing DataFrames
    train_df = train_df[train_df['pft'] != "other_herb"]
    train_df = train_df.drop(["shapeID"],axis=1) #'site', 'year'
    test_df = test_df[test_df['pft'] != "other_herb"]
    test_df = test_df.drop(["shapeID"],axis=1) #'site','year'
    # check the class distribution in the training and testing sets
    print("Training set class distribution:")
    print(train_df['pft'].value_counts(normalize=True))
    print("Testing set class distribution:")
    print(test_df['pft'].value_counts(normalize=True))

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

    rf_model_path = os.path.join(rf_output_dir, f"rf_model_{shapefile_description}.joblib")
    if not os.path.exists(rf_model_path):
        rf_model = RF(n_estimators=ntree, random_state=42, class_weight="balanced")
        rf_model.fit(train_df.drop("pft", axis=1), train_df["pft"])
        dump(rf_model, rf_model_path)
    else:
        rf_model = load(rf_model_path)

    y_pred_val = rf_model.predict(test_df.drop("pft", axis=1))
    accuracy_val = accuracy_score(test_df["pft"], y_pred_val)
    print("Validation Accuracy:", accuracy_val)
    print("Validation Classification Report:")
    print(classification_report(test_df["pft"], y_pred_val))
    print("Validation Confusion Matrix:")
    print(confusion_matrix(test_df["pft"], y_pred_val))
    
    # # Visualize classifications
    # # in the testing dataframe, group to identify which plants are near each other
    # # create a new column called "groupID"
    # vis_predict_df = pd.concat([test_df.drop(['index'],axis=1).reset_index(drop=True),pd.DataFrame(y_pred_val,columns=['pft_pred'])],axis=1,join='outer')
    # vis_predict_df['groupID'] = np.nan
    # group_id = 0
    # # loop through each row in the DataFrame
    # for i, row in vis_predict_df.iterrows():
    #     # check if the row has already been assigned a group ID
    #     if np.isnan(row['groupID']):
    #         group_id += 1
    #         vis_predict_df.loc[i, 'groupID'] = group_id

    #         # find all rows within 30m of the current row
    #         nearby_rows = vis_predict_df[(vis_predict_df['eastingIDs'] >= row['eastingIDs'] - 30) &
    #                         (vis_predict_df['eastingIDs'] <= row['eastingIDs'] + 30) &
    #                         (vis_predict_df['northingIDs'] >= row['northingIDs'] - 30) &
    #                         (vis_predict_df['northingIDs'] <= row['northingIDs'] + 30)]
    #         # assign the nearby rows to the same group ID
    #         vis_predict_df.loc[nearby_rows.index, 'groupID'] = group_id

    # # Create a colormap with a discrete color ramp
    # cmap = cm.get_cmap('tab10', len(vis_predict_df['pft'].unique()))
    # color_dict = {value: cmap(i) for i, value in enumerate(vis_predict_df['pft'].unique())}    
    # vis_predict_df['color'] = vis_predict_df['pft'].map(color_dict) # Add a new column with the corresponding colors

    # # Plot predictions
    # for i in vis_predict_df['groupID'].unique():
    #     vis_predict_temp = vis_predict_df[vis_predict_df['groupID']==i]
    #     fig, ax = plt.subplots(figsize=(10, 10))

    #     x = int(vis_predict_temp['eastingIDs'].unique())
    #     y = int(vis_predict_temp['northingIDs'].unique())
        
    #     xmin = x+1000
    #     xmax = x
    #     ymin = y
    #     ymax = y-1000
    #     for i, row in vis_predict_temp.iterrows():
    #         if (x+serial_index_to_coordinates(row['pixelNumber'],1000)[0]) < xmin:
    #             xmin = x+serial_index_to_coordinates(row['pixelNumber'],1000)[0]
    #         if (x+serial_index_to_coordinates(row['pixelNumber'],1000)[0]) > xmax:
    #             xmax = x+serial_index_to_coordinates(row['pixelNumber'],1000)[0]
    #         if (y-(1000-serial_index_to_coordinates(row['pixelNumber'],1000)[1])) < ymin:
    #             ymin = y-(1000-serial_index_to_coordinates(row['pixelNumber'],1000)[1])
    #         if (y-(1000-serial_index_to_coordinates(row['pixelNumber'],1000)[1])) > ymax:
    #             ymax = y-(1000-serial_index_to_coordinates(row['pixelNumber'],1000)[1])

    #     # search for the image that includes the specified easting and northing coordinate
    #     imagery_folder = os.path.join(data_int_path,row['site'],str(row['year']),"stacked_aop")
    #     for file in os.listdir(imagery_folder):
    #         if file.endswith('.tif'):
    #             imagery_path = os.path.join(imagery_folder, file)
    #             imagery = rasterio.open(imagery_path)
    #             if imagery.bounds.left == x and y == imagery.bounds.top:

    #                 # plot the remote sensing pixels used for validation with a 10 pixel buffer
    #                 window = rasterio.windows.Window(xmin-10, ymin-10, xmax-xmin+20,ymax-ymin+20)
    #                 data = imagery.read(window=window)
    #                 # multiband image
    #                 # specify which three bands to plot
    #                 bands = [9, 10, 11]  # e.g., plot bands 4, 3, and 2 as RGB                   
    #                 data = imagery.read(bands)  # read the specified bands                    
    #                 data = np.transpose(data, (1, 2, 0)) # transpose the data to the correct shape for plotting

    #                 # plot the data with RGB colors
    #                 ax.imshow(data, cmap='viridis', extent=(xmin-10, xmax+10, ymin-10, ymax+10))

    #                 # loop through each row in the DataFrame
    #                 for i, row in vis_predict_temp.iterrows():
    #                     x_temp = serial_index_to_coordinates(row['pixelNumber'],1000)[0]
    #                     y_temp = serial_index_to_coordinates(row['pixelNumber'],1000)[1]
                        
    #                     # plot the classified predictions on top with a small amount of transparency
    #                     ax.scatter(x, y, c=row['color'], alpha=0.5, label=row['pft']) 

    #                 break
         
    #     ax.legend() # add a legend
    #     # save the plot
    #     plt.savefig(os.path.join(rf_output_dir, "validation", row['site']+"_"+str(row['year'])+'_'+str(round(imagery.bounds.left))+'_'+str(round(imagery.bounds.top))+'.png'))
    
    log.info('Trained PFT classifier saved in this folder: '
             f'{rf_model_path}')
    
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


