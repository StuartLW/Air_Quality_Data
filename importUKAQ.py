import rdata
from urllib.request import urlopen
import pandas as pd
import numpy as np
from tqdm import tqdm
import itertools
import h5py

# constants from open air project
META_URLS= {'aurn':"https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData",
            'saqn' : "https://www.scottishairquality.scot/openair/R_data/SCOT_metadata.RData",
            'saqd' : "https://www.scottishairquality.scot/openair/R_data/SCOT_metadata.RData",
            'ni' : "https://www.airqualityni.co.uk/openair/R_data/NI_metadata.RData",
            'waqn' : "https://airquality.gov.wales/sites/default/files/openair/R_data/WAQ_metadata.RData",
            'aqe' : "https://airqualityengland.co.uk/assets/openair/R_data/AQE_metadata.RData",
            'local' : "https://uk-air.defra.gov.uk/openair/LMAM/R_data/LMAM_metadata.RData",
            'lmam' : "https://uk-air.defra.gov.uk/openair/LMAM/R_data/LMAM_metadata.RData"
            }
DATA_URL_ROOTS= {'aurn': 'https://uk-air.defra.gov.uk/openair/R_data/',
                 "aqe" : "https://airqualityengland.co.uk/assets/openair/R_data/",
                "saqn" : "https://www.scottishairquality.scot/openair/R_data/",
                "waqn" : "https://airquality.gov.wales/sites/default/files/openair/R_data/",
                "ni" : "https://www.airqualityni.co.uk/openair/R_data/",
                "local" : "https://uk-air.defra.gov.uk/openair/LMAM/R_data/"
                }
DATA_TYPES= {'hourly':'',
             '15min':'_15min',
             'daily':'_daily_mean',
             '8_hour':'_8hour_mean',
             '24_hour':'_24hour_mean',
             'daily_max_8':'_daily_max_8hour'
             }
# dictionary to hold meta data
meta_data = {}

def get_metadata(source='aurn', force_update = False)->None:
    """get the meta data file for use in the project. Meta data files contain information on the location of different sites as well
    as the data available and the dates that data acquisition has taken place between.  This should not need to be fetched more than
    once in a session therefore the data is only fetched if it does not exist, unless force_update is set to True.
    end date of 'ongoing is set to current date and an adidtional ongoing flag column set
    ratified_to values of 'Never' are set to 1-1-1900 to enable date columns to be handled as datetime objects

    Args:
        source (str, optional): The network required. Defaults to 'aurn'.
        force_update (bool, optional): Whether to force an update of the network metadata. Defaults to False.

    """
    if source not in meta_data.keys()  or force_update:
        if source in META_URLS.keys():
            meta_dict = _get_RData(META_URLS[source])
            meta_data[source]= meta_dict[list(meta_dict.keys())[0]]
            # clean up dates and ongoing column
            meta_data[source]['ongoing_flag']= meta_data[source]['end_date'].map(lambda x: True if x=='ongoing' else False)
            # set end date to current date if ongoing
            meta_data[source]['end_date']= meta_data[source]['end_date'].map(lambda x: pd.Timestamp.today().strftime('%Y-%m-%d') if x=='ongoing' else x)
            # set ratified to column to 1-1-1900 if Never
            meta_data[source]['ratified_to']= meta_data[source]['ratified_to'].map(lambda x: '1900-01-01' if x=='Never' else x)
            # convert date columns to datetime objects
            meta_data[source][['start_date', 'end_date', 'ratified_to']]=meta_data[source][['start_date', 'end_date', 'ratified_to']].apply(pd.to_datetime)

def get_sites(source)->pd.DataFrame:
    """network metadata contains multiple entries for each site as it lists the pollutants measured. This function returns a dataframe the location, and id information only.
    Note that the presence of a site in the list does nto mean data will be available as the years of operation and the pollutants measuredd at each site differ.

    Args:
        source (str, optional): The source network for the data.

    Returns:
        pd.DataFrame: A data frame of the core locationa ttributes of the sites
    """
    if source not in meta_data.keys():
        get_metadata(source)
    df=meta_data[source][['site_id', 'site_name', 'location_type', 'latitude', 'longitude','zone', 'agglomeration', 'local_authority']].drop_duplicates()
    return df

def get_parameters(source)->list:
    """returns a list of available parameter names. Note not all parameters will be available at all sites

    Args:
        source (str): the source to query

    Returns:
        list: a list of avaiable parameter names
    """
    if source not in meta_data:
        get_metadata(source)
    parameters = list(meta_data[source]['parameter'].unique())
    return parameters

def guess_all(year, source, pollutant_list):
    """for a given source, year range and pollutant list returns the list of sites with data for the given range.  For readign hdf file on nly the year is checked not the pollutants

    Args:
        year (list): years to be covered
        source (str): data source identifier eg 'aurn'
        pollutant_list (list): list of target polutants

    Returns:
        list: list of available sites with required data
    """

    if source in META_URLS:
        # get the metadata if we have not already done it
        if source not in meta_data:
            get_metadata(source)
        md=meta_data[source]
        # find the sites with any data acquired over the period
        filtered_ids = set(md[(md['start_date'].dt.year<=min(year)) &(md['end_date'].dt.year>=max(year))]['site_id'])
        if 'any' not in pollutant_list:
            # add additional filters based on polutants
            for  pollutant in pollutant_list:
                # check if data is available for that pollutant in the time window
                pollutant_site_ids = set(md[(md['start_date'].dt.year<=min(year)) &(md['end_date'].dt.year>=max(year)) & (md['parameter'] == pollutant)]['site_id'])
                # find interesection
                filtered_ids = filtered_ids & pollutant_site_ids
    else: # hdf file
        print("guesses of 'all' for sites in an hdf file are limited to years where data is present, it is currently not possible to check individual pollutants")
        hdf_file = h5py.File(source)
        keys = hdf_file.keys()
        sources, years = zip(*[[x.split('_')[0],int(x.split('_')[1])] for x in list(keys)])
        sources= np.array(sources)
        years= np.array(years)
        filter = np.logical_and(years>=min(year),years<=max(year))
        filtered_ids = list(map(str,np.unique(sources[filter])))
        hdf_file.close()
    return list(filtered_ids)

def importUKAQ(site,
                year,
                source, 
                data_type='hourly',
                pollutant = 'any',
                to_narrow= False,
                verbose=False,
                progress=True,
                write_raw_toHDF= False
                )->pd.DataFrame:
    """Downloads and parses RData files from open air quaity networks into a dataframe.  If multiple monitoring sites are specified along
    with multiple years all years for all sites are downloaded. Also contains options to save raw downloads to an hdf5 format which can
    be reloaded with teh same function.  Which subtype of data is to be loaded can also be defined

    Args:
        site (str or list): A code or list of codes representing the monitoring site(s) eg my1. Alternatively 'all' will get all sites that meet the year and pollutant criteria
        If more than one source is being used then the site list will attempt to be coerced in a manner that the sites are linked to teh source.site 1 and site 2 from source 1 
        and all sites from source2
        year (int, str or list): The year(s) to be recovered. Can be a single year, a list of years. Alternatively
        a string of eg '2000:2005' will be translated into an inclusive list.
        source (str or list): the source of the data, eg 'aurn'.  Alternatively set to a local file path to read from local hdf5 file, previously created.
        If a list of sources is provided then the site argument must also be a list of the same length (which can be a list of lists if more than one site is needed from each source)
        data_type (str, optional): The data type (frequency of observation) to recover. Not all frequencies are available from all sources. Defaults to 'hourly'.
        pollutant (str or list, optional) - list of pollutants to . defaults to 'any'be returned. Also used in guessign sites if 'all' sites are specified.  A value of 'any' keeps all data. Note all data is downloaded regardless
        of this value - the selection takes place after download. Defaults to 'any;
        to_narrow (bool, optional). Melt the output dataframe so that the dataframe is returned with a
        column identifying the pollutant name and a column containing the corresponding concentration/statistic. Defaults to False
        verbose (bool, optional): if True will print which data is being downloaded Defaults to False.
        progress (bool, optional): Shows a progress bar for downloads. Defaults to True.
        write_raw_toHDF (bool or string, optional): If a string file path is provided then an hdf5 file containing all the downloaded data will
        be written locally.Useful to enable rereading of all data without redownloading Defaults to False.

    Raises:
        TypeError: for incorrectly defined sites or year
        ValueError: if source or data type not allowed values

    Returns:
        pd.DataFrame: concatenated pandas dataframe of the data_type requested
    """
    # parse the arguments into a common format
 
    # year can be an integer, list of intgers or a string range
    if isinstance(year, int):
        year_list = [year]
    elif isinstance(year,str):
        # split a range string and create a list of years
        year_list = list(range(int(year.split(':')[0]),int(year.split(':')[1])+1))
    elif isinstance(year, list):
        year_list = year
    else:
        raise TypeError("year must be an integer, a list of integers or a string range wih format eg '2005:2010'")
    

    # process pollutants
    if isinstance(pollutant,str):
        pollutant_list = [pollutant]
    else:
        pollutant_list = pollutant


    # process the source_list we need each source to be a list for the matching with potentially many sites
    if isinstance(source,str):
        source_list = [source]
    elif isinstance(source, list):
        source_list = source
    # parse source to url , unless local file
    source_url_list = [[DATA_URL_ROOTS[x]]  if x in DATA_URL_ROOTS.keys() else [x] for x in source_list]
    source_list = [[x] for x in source_list]
   
    # Parse the sites. Need to handle the possibility of sites being listed as 'all' and also there to be more than one source
    # A bit clunky and could probably have better error checking    
    if len(source_list) ==1: # only a singe site so simple treatment
            if isinstance(site,str):
            # single site
                site_list = [_parse_sites(site, year_list, source_list[0], pollutant_list, var_type='str')]
            else:
                site_list = [[_parse_sites(s,year_list, source_list[0],pollutant_list) for s in site]]
    elif len(source_list) == len(site):
        #  we have multiple sites and multiple entries for sources as required, need to handle mixture of 'all' all and site name
        # 'all cannot appear in second level list so should not be parsed in tha else case
        site_list = []
        for x,y in zip(site, source_list):
            if isinstance(x,str):
                site_list += [_parse_sites(x, year_list, y, pollutant_list,var_type='str')]
            else:
                site_list.append( [z.upper() for z in x])
        # site_list = [_parse_sites(x, year_list, y[0], pollutant_list,var_type='str') if isinstance(x,str) else [z.upper() for z in x] for x,y in zip(site, source_list)]
    else:
        raise ValueError(f'There are {len(source_list)} sources specified but the site description cannot be coerced in manner that enable correct linkage of site and source')


    if data_type not in DATA_TYPES.keys():
        raise ValueError('data_type must be one of {}'.format(DATA_TYPES.keys())) 
    
    
    source_site_list = [list(itertools.product(x,y)) for x,y in zip(source_url_list,site_list,)]
    source_site_list = [element for innerList in source_site_list for element in innerList]
    # create all the permutations of site and year
    job_list = list(itertools.product(source_site_list,year_list))

   
    results={} # an empty dictionary to hold the results

    # main loop to get the data file
    for job_params in tqdm(job_list, disable=not progress):
        if job_params[0][0][:4]=='http':
            url = f'{job_params[0][0]}{job_params[0][1]}_{job_params[1]}.RData'
            if verbose:
                print(f'getting data for {job_params[0][1]} for {job_params[1]}')
            try:
                results[f'{job_params[0][1]}_{job_params[1]}']= _get_RData(url)
            except ValueError:
                print("{url} url does not exist, probably because data does not exist for this site and year combination")
        else: # read from local hdf file
            try:
                results[f'{job_params[0][1]}_{job_params[1]}']= _get_HDFData(job_params)
            except ValueError:
                print(f"{job_params}: data does not exist, probably because data does not exist for this site and year combination")
             
   # write raw files to HDF if file path specified - to enable rereading faster
    if write_raw_toHDF != False:
        _write_raw_files_toHDF(results, write_raw_toHDF)
    
    # process the data to create a dataframe
    selected_results = [results[k][f'{k}{DATA_TYPES[data_type]}'] for k in results.keys()]
    out_df=pd.concat(selected_results, axis=0)
    # fix formatting and onvert to date time object
    out_df['site'] = out_df['site'].astype('str')
    out_df['date'] = pd.to_datetime(out_df['date'], unit='s')
    if 'any' not in pollutant_list:
        # only keep data from the pollutant list for the numeric data types
        cols_to_drop = np.setdiff1d(out_df.select_dtypes(include='number').columns, pollutant_list)
        out_df.drop(cols_to_drop, axis=1, inplace = True)
    if to_narrow:
        out_df = out_df.melt(id_vars=['site','code','date'])
    # clean up
    del results
       
    return out_df

def _parse_sites(site, year,source,pollutant, var_type='list_member'):
    """ helper function to parse site arguments

    Args:
        site (str|list): the site parameter
        year (list): the generated year list
        source (list): the generated list of sites
        pollutant (str|list): the list of desired polutants
        var_type (str, optional): a flag for whetehr we are parsing a string or a list.  Defaults to 'list_member'.

    Returns:
        list: a list of sites in the correct format
    """

    if site == 'all':
         rtn = guess_all(year, source[0], pollutant)
    else:
        if var_type=='str':
            rtn = [site.upper()]
        else:
            rtn = site.upper()
                         
    return rtn


def _get_RData(url)->dict:
    """gets an Rdata file from an open air source and returns a dictionary of dataframes

    Args:
        url (str): The url to get the Data from - 

    Returns:
        dict: A dctionary object containing dataframes
    """
    try:
        with urlopen(url) as dataset:
            parsed = rdata.parser.parse_file(dataset)
        converted = rdata.conversion.convert(parsed)
    except:
        raise ValueError("Cant open {}, probably does not exist".format(url))
    return converted

def _get_HDFData(job_param)->dict:
    """ Gets data from an hdf file of raw data previously created from downloads to give the same 
    data format as if the files had been downloaded again - enables local storage

    Args:
        job_param (tuple): tuple of ((file, site_id), year)

    Raises:
        ValueError: if data defined by the tuple does not exist

    Returns:
        dict: dictionary of dataframes of different data intervals - equivalent to that generated from download
    """
    prefix = f'{job_param[0][1]}_{job_param[1]}'
    hdf_file = h5py.File(job_param[0][0])
    if prefix in hdf_file.keys():
        res_dict = {}
        for key in hdf_file[prefix].keys():
            res_dict[key]= pd.read_hdf(job_param[0][0], key=f'{prefix}/{key}')
    else:
        raise ValueError("specified data not in hdf file")
    hdf_file.close()
    return res_dict
   
def _write_raw_files_toHDF(results, file)->None:
    """writes the data from the raw download to an hdf5 file with the same structure as the original downloaded RData

    Args:
        results (dict): dictionary of downloaded data frames
        file (str): filepath to the hdf5 file 
    """
    for key in results.keys():
        for data_set in results[key].keys():
            df = results[key][data_set]
            # need to convert site to string to enable saving
            df['site'] =df['site'].astype('str')
            df.to_hdf(file, key = f'{key}/{data_set}',format='table', mode='a')