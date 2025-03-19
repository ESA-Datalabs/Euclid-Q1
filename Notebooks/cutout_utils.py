from astroquery.esa.euclid.core import EuclidClass
from tqdm._tqdm_notebook import tqdm_notebook
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D
from skimage.transform import resize
from astropy.io import fits
import astropy.units as u
from astropy import wcs
import numpy as np
import re

#choosing environment
Euclid = EuclidClass(environment='PDR')


###################################### Functions for matching the sources with files ######################################
            
def chooseBestFile(df_res_dup, radec_colnames):
    """Choosing best file for sources with multiple by checking where the source is closest to image centre"""
    
    min_dist_per_plotcode = df_res_dup.groupby([radec_colnames["ra_colname"], 'filter_name'], as_index=False, sort=False)['dist'].idxmin()
    df_final = df_res_dup.loc[min_dist_per_plotcode["dist"]]
    return df_final


def addFilters(product_type, query, instrument_name, nisp_filters=[]):
    """In case nisp_filters is specified modify query to get files for only those filters"""
    
    filter_part = ""
    query_parts = query.split("ON ")
    for one_filter in nisp_filters:
        if one_filter != nisp_filters[-1]:
            filter_part += "("+ product_type + ".filter_name='" + one_filter + "') OR "
        else:
            filter_part += "(" + product_type + ".filter_name='" + one_filter + "')"
    query = query_parts[0] + "ON (" + filter_part + ") AND" + query_parts[1]
    return query

def getRaDecColnames(all_colnames):
    """Return column names used for ra, dec"""
    
    if 'ra' in all_colnames and 'dec' in all_colnames:
        ra_colname, dec_colname = 'ra', 'dec'
    else:
        ra_colname, dec_colname = 'right_ascension', 'declination'
    return {"ra_colname": ra_colname, "dec_colname": dec_colname}


def checkColnames(source_query):
    """Load metadata for user table and get column names for ra, dec by calling getRaDecColnames"""
    
    table = Euclid.load_table(table=source_query, verbose=False)
    all_colnames = [col.name for col in table.columns]
    return getRaDecColnames(all_colnames)


def getFilesCalib(source_query, instrument_name, nisp_filters, searchdist, radec_colnames):
    """Match sources with calibrated files by checking if file fov contains source"""
    
    query= f"""SELECT sources.*, frames.file_name, frames.file_path, frames.datalabs_path, frames.instrument_name, frames.filter_name, frames.observation_id, frames.ra AS image_ra, frames.dec AS image_dec, DISTANCE(frames.ra, frames.dec, sources.{radec_colnames["ra_colname"]}, sources.{radec_colnames["dec_colname"]}) AS dist
    FROM {source_query} as sources
    JOIN sedm.calibrated_frame AS frames
    ON (frames.instrument_name='{instrument_name}') AND (frames.fov IS NOT NULL AND CONTAINS(CIRCLE('ICRS', sources.{radec_colnames["ra_colname"]}, sources.{radec_colnames["dec_colname"]}, {searchdist/60}), frames.fov)=1)"""

    #modify query to include only the provided nisp filters
    if instrument_name == 'NISP' and len(nisp_filters) > 0 and len(nisp_filters) < 3:
        query = addFilters("frames", query, instrument_name, nisp_filters)
        
    #make query to get file info
    job = Euclid.launch_job_async(query, verbose=False)
    
    if job == None:
        raise AttributeError("Query for files failed. Please make sure your source query or user table contains columns named 'right_ascension' and 'declination' or 'ra' and 'dec'. When using a source query it should have brackets around it.")
        
    df_res = job.get_results().to_pandas()
    print("Found", len(df_res), "query results")
    
    return df_res
    
def getFilesMosaic(source_query, instrument_name, nisp_filters, radec_colnames):
    """Match sources with mosaic files by checking if file tile_index matches the first 9 digits of source's segmentation_map_id"""
    
    query= f"""SELECT sources.*, mosaics.file_name, mosaics.file_path, mosaics.datalabs_path, mosaics.mosaic_product_oid, mosaics.tile_index, mosaics.instrument_name, mosaics.filter_name, mosaics.ra AS image_ra, mosaics.dec AS image_dec, DISTANCE(mosaics.ra, mosaics.dec, sources.{radec_colnames["ra_colname"]}, sources.{radec_colnames["dec_colname"]}) AS dist
        FROM {source_query} as sources
        JOIN sedm.mosaic_product AS mosaics
        ON (mosaics.instrument_name='{instrument_name}') AND (CAST(mosaics.tile_index AS CHAR(25)) = SUBSTRING(CAST(segmentation_map_id AS CHAR(25)), 1, 9))"""

    #modify query to include only the provided nisp filters
    if instrument_name == 'NISP' and len(nisp_filters) > 0 and len(nisp_filters) < 3:
        query = addFilters("mosaics", query, instrument_name, nisp_filters)
        
    #make query to get file info
    job = Euclid.launch_job_async(query, verbose=False)
    
    if job == None:
        raise AttributeError("Query for files failed. Please make sure your source query or user table contains columns named 'right_ascension' and 'declination' or 'ra' and 'dec'. When using a source query it should have brackets around it.")
    
    df_res = job.get_results().to_pandas()
    print("Found", len(df_res), "query results")
    
    return df_res


def getFilesMosaicBackup(source_query, instrument_name, nisp_filters, searchdist, radec_colnames):
    """Match sources with mosaic files by checking if file fov contains source"""
    
    query= f"""SELECT sources.*, mosaics.file_name, mosaics.file_path, mosaics.datalabs_path, mosaics.mosaic_product_oid, mosaics.tile_index, mosaics.instrument_name, mosaics.filter_name, mosaics.ra AS image_ra, mosaics.dec AS image_dec, DISTANCE(mosaics.ra, mosaics.dec, sources.{radec_colnames["ra_colname"]}, sources.{radec_colnames["dec_colname"]}) AS dist
        FROM {source_query} as sources
        JOIN sedm.mosaic_product AS mosaics
        ON (mosaics.instrument_name='{instrument_name}') AND (mosaics.fov IS NOT NULL AND CONTAINS(CIRCLE('ICRS', sources.{radec_colnames["ra_colname"]}, sources.{radec_colnames["dec_colname"]}, {searchdist/60}), mosaics.fov)=1)"""

    #modify query to include only the provided nisp filters
    if instrument_name == 'NISP' and len(nisp_filters) > 0 and len(nisp_filters) < 3:
        query = addFilters("mosaics", query, instrument_name, nisp_filters)
    
    #make query to get file info
    job = Euclid.launch_job_async(query, verbose=False)
    
    if job == None:
        raise AttributeError("Query for files failed. Please make sure your source query or user table contains columns named 'right_ascension' and 'declination' or 'ra' and 'dec'. When using a source query it should have brackets around it.")
    
    df_res_dup = job.get_results().to_pandas()
    df_res = chooseBestFile(df_res_dup, radec_colnames)
    print("Found", len(df_res), "query results")
    
    return df_res

def checkInput(instrument_name, nisp_filters):
    """Checking user input value correctness for instrument name and nisp filters"""
    
    valid_instruments = ["VIS", "NISP"]
    valid_filters = ["NIR_H", "NIR_J", "NIR_Y", ""]
    if instrument_name not in valid_instruments:
        raise ValueError(f"Invalid instrument name. Choose from: {valid_instruments}") 
    for filter_name in nisp_filters:
        if filter_name not in valid_filters:
            raise ValueError(f"Invalid filter name. Choose from: {valid_filters}")     
    return

def getFiles(product_type, sources, instrument_name, nisp_filters=[], segmentation_map_id = False, searchdist = 0.5):
    """Matching sources with files to make cutouts from
    
    Args:
        product_type (str, mandatory): either mosaic or calibrated_frame
        sources (str, mandatory): either the query constructed above as str or the name of the uploaded user table
        instrument_name (str, mandatory): either 'VIS' or 'NISP'(specifying NISP will get files for all 3 bands unless nisp_filters specified)
        nisp_filters (list, optional): list of strings containing any or all of the NISP filter names NIR_H, NIR_J, NIR_Y
        segmentation_map_id (bool, optional): whether the source table has this column (for faster file matching for mosaics, default False)
    """
    
    #checking input value correctness
    checkInput(instrument_name, nisp_filters)
    
    #archive login for queries
    Euclid.login(credentials_file = '/media/user/cred.txt')
    
    # check if user table uses columns ra, dec or right_ascension, declination
    if "user_" in sources:
        radec_colnames = checkColnames(sources)
    else:
        radec_colnames = {"ra_colname": "right_ascension", "dec_colname": "declination"}
    
        
    #getting the files for sources
    if product_type == "mosaic":
        if segmentation_map_id == True:
            df = getFilesMosaic(sources, instrument_name, nisp_filters, radec_colnames)
        else:
            df = getFilesMosaicBackup(sources, instrument_name, nisp_filters, searchdist, radec_colnames)        
    elif product_type == "calibrated_frame":
        df = getFilesCalib(sources, instrument_name, nisp_filters, searchdist, radec_colnames)
    else:
        raise ValueError(f"Invalid instrument name. Choose from: mosaic, calibrated_frame") 
        
    return df

###################################### Functions for making the cutouts ######################################

def getDetIdx(hdul, ra, dec, mean_gap_x, mean_gap_y, det_width, det_height):
    """Calculating the index of the detector the source is on"""
    
    #getting source pixel coordinates in the frame of the first detector
    w1 = wcs.WCS(hdul[1].header)
    px, py = w1.wcs_world2pix(ra, dec, 1)
    
    #if negative in the frame of the first detector something is wrong
    if px < 0 or py < 0:
        return -1, -1
    
    #moving reference frame so the transitions are in the middle of the gaps
    px = px + mean_gap_x / 2
    py = py + mean_gap_y / 2
    
    #calculate detector index knowing gap and detector size
    det_idx_x = int(np.ceil(px / (det_width + mean_gap_x)))
    det_idx_y = int(np.ceil(py / (det_height + mean_gap_y)))
    
    return det_idx_x, det_idx_y

def getExtIdxVIS(hdul, ra, dec):
    """Getting the index of the fits extension the source is on for VIS image"""
    
    #average gaps between detectors in either direction
    mean_gap_x, mean_gap_y = 134.7317209303794, 612.8441331180677
    det_x, det_y = 2048, 2066
    
    #getting the detector index (11 to 66)
    det_idx_x, det_idx_y = getDetIdx(hdul, ra, dec, mean_gap_x, mean_gap_y, 2*det_x, 2*det_y)
    
    if det_idx_x > 6 or det_idx_x < 0 or det_idx_y > 6 or det_idx_y < 0:
        return -1
    
    #getting the fits extension number of the E quadrant of the detector from the index of the detector
    ext_idx_E = hdul.index_of(str(det_idx_y) + "-" + str(det_idx_x) + ".E.SCI")
    
    #finding the correct quadrant of the detector (checking all 4)
    for ext in range(ext_idx_E, ext_idx_E+12, 3):
        w = wcs.WCS(hdul[ext].header)
        px, py = w.wcs_world2pix(ra, dec, 1)
        if py < det_y and py > 0 and px < det_x and px > 0:
            return ext
        
    return -1

def getExtIdxNISP(hdul, ra, dec):
    """Getting the index of the fits extension the source is on for NISP image"""
    
    #average gaps between detectors in either direction
    mean_gap_x, mean_gap_y = 164.90065775050047, 328.63495367441845
    det_dim = 2040
    
    #getting the detector index (11 to 44)
    det_idx_x, det_idx_y = getDetIdx(hdul, ra, dec, mean_gap_x, mean_gap_y, 2040, 2040)
    
    if det_idx_x > 4 or det_idx_x < 0 or det_idx_y > 4 or det_idx_y < 0:
        return -1
    
    #getting the fits extension number from the index of the detector
    ext_idx = hdul.index_of("DET" + str(det_idx_y) + str(det_idx_x) + ".SCI")
    
    #checking that source is not in gap
    w = wcs.WCS(hdul[ext_idx].header)
    px, py = w.wcs_world2pix(ra, dec, 1)
    if py < det_dim and py > 0 and px < det_dim and px > 0:
        return ext_idx
    
    return -1

def cutoutsGroup(ra, dec, file_name, hdul, product_type, cutout_size_inf, cutout_func, optional_cutout_params, output_location, resize_params):
    """Making a cutout for a source"""
    
    # if calibrated file then check which fits file extension to cut otherwise for mosaics cut primery extension in fits file
    if product_type == "frame":
        if "_NIR_" in file_name:
            ext_nr = getExtIdxNISP(hdul, ra, dec)
        elif "_VIS_" in file_name:
            ext_nr = getExtIdxVIS(hdul, ra, dec)
        else:
            ext_nr = -1
            

        if ext_nr == -1:
            # if no detector found (source in gap) then don't try to cut
            return 0
    else:
        ext_nr = 0
        
    # making cutout
    data = hdul[ext_nr].data
    header = hdul[ext_nr].header
    
    if cutout_func != None:
        cutout = Cutout2D(data= data, position = SkyCoord(ra, dec, frame='icrs', unit="deg"), size=cutout_func(cutout_size_inf), wcs=wcs.WCS(header), **optional_cutout_params)
    else:
        cutout = Cutout2D(data=data, position = SkyCoord(ra, dec, frame='icrs', unit="deg"), size= cutout_size_inf, wcs=wcs.WCS(header), **optional_cutout_params)
    
    cutout_data = cutout.data
    
    #resizing
    if resize_params != None:
        cutout_data = resize(cutout_data, **resize_params)
    
    hdu_cutout = fits.PrimaryHDU(cutout_data)
    hdu_cutout.header.update(cutout.wcs.to_header())
    
    #naming and saving
    if product_type == "mosaic":
        name = re.search("\-(.*)\-", file_name).group()[1:-1] + "-CUTOUT_" + str(round(ra, 7)) + "_" + str(round(dec, 7)) + ".fits"
    elif product_type == "frame":
        name = re.search("\_(.*)\_", file_name).group()[1:-1] + "-CUTOUT_" + str(round(ra, 7)) + "_" + str(round(dec, 7)) + ".fits"
        
    hdu_cutout.writeto(output_location + name, overwrite=True)
    print("Made cutout", name)
    
    return 1

def cuotutSubdf(sub_df, product_type, cutout_size, cutout_func, optional_cutout_params, output_location, resize_params, radec_colnames):
    """Read in image once and make cutouts for each row/source in subdataframe by calling cutoutsGroup function"""
    
    #load image once
    nr_of_cutouts = 0
    
    # getting the image from volume
    hdul = fits.open(sub_df["datalabs_path"].iloc[0] + '/'+ sub_df["file_name"].iloc[0])
    
    #check if doing constant radius cutouts or variable and apply cutout function to each row in small df
    if type(cutout_size) == str:
        ap_res = sub_df.apply(lambda x: cutoutsGroup(x[radec_colnames["ra_colname"]], x[radec_colnames["dec_colname"]], x["file_name"], hdul, product_type, x[cutout_size], cutout_func, optional_cutout_params, output_location, resize_params), axis=1)
    else:
        ap_res = sub_df.apply(lambda x: cutoutsGroup(x[radec_colnames["ra_colname"]], x[radec_colnames["dec_colname"]], x["file_name"], hdul, product_type, cutout_size, cutout_func, optional_cutout_params, output_location, resize_params), axis=1)
        
    nr_of_cutouts += ap_res.sum()
    return nr_of_cutouts
    
def makeCutoutsByTile(df, output_location, cutout_size, cutout_func=None, optional_cutout_params = {}, resize_params=None):
    """Check inputs and make cutouts for the provided df by first grouping the df by file name and calling cuotutSubdf function for each group

    Args:
        df (astropy table, mandatory): table that includes sources and their file info - output of getFiles function
        output_location (str, mandatory): path where to store the cutouts
        cutout_size (int, tuple, str, mandatory): input size to astropy.Cutout2D. Can be anything astropy Cutout2D allows or in case of variable size the name of the column from the source table that the cutout size should depend on.
        cutout_func (mandatory if cutout_size is str)- in case of variable cutout size this is the name of a user defined function how the size depends on the selected column. Default None.
        optional_cutout_params (dict, optional): other inputs to astropy.Cutout2D like mode and fill value. The default mode is trim so when a source is on the edge of the image the cutout will be smaller.
        resize_params (dict, optional) - in case you want to resize the cutout before saving it then the inputs to skimage resize. By default no resizing is done.

    """
    
    # check whether column names are right_ascension + declination or ra + dec
    radec_colnames = getRaDecColnames(list(df.columns))
    
    
    #checking if doing mosaic cutouts or calibrated frame ones from the first file name in df
    if "MOSAIC" in df["file_name"].iloc[0]:
        product_type = 'mosaic' 
    else:
        product_type = 'frame'
    
    #if inputs are correct group the df by file name and apply cutout funtion to each subdf
    if (isinstance(cutout_size, str) and cutout_func!=None) or ((isinstance(cutout_size, u.quantity.Quantity) or isinstance(cutout_size, int) or isinstance(cutout_size, tuple)) and cutout_func==None):
        if ((resize_params == None) or (resize_params != None and "output_shape" in resize_params)):
            res = [cuotutSubdf(sub_df, product_type, cutout_size, cutout_func, optional_cutout_params, output_location, resize_params, radec_colnames) for index, sub_df in tqdm_notebook(df.groupby(['file_name']), desc = "Number of files processed")]
            print("Created a total of", sum(res), "cutouts")
        else:
            raise ValueError(f"Invalid set of input input arguments - resize_params (if given) must include output_shape.")
    else: 
        raise ValueError(f"Invalid set of input input arguments for cutout radius. If you are trying to use a constant radius set cutout_radius as a astropy.units.quantity.Quantity type object and don't specify the cutout_func. If you are trying to use a variable radius then specify cutout_radius as a string corresponding to a column name in the dataframe and also provide a function name to cutout_func to specify how to make the radious dependent on that source characteristic.")
    return 