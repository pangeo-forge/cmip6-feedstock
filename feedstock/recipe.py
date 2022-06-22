# Refactored logic
# Key assumptions: 
# - We ever only request a single IID. Broader searches are implemented in the parsing module (To be commited)
# - We start with an instance id, but will actually search for a dataset id with 'latest' query set to true. This means if an outdated iid is requested we will error out (TODO: Need to think about how to handle this)

#TODO: Check multiple nodes with some sort of priority (naomi had done quite a bit with that)
# node_pref = {
#     "esgf-data1.llnl.gov": 0,
#     "esgf-data2.llnl.gov": 0,
#     "aims3.llnl.gov": 0,
#     "esgdata.gfdl.noaa.gov": 10,
#     "esgf-data.ucar.edu": 10,
#     "dpesgf03.nccs.nasa.gov": 5,
#     "crd-esgf-drc.ec.gc.ca": 6,
#     "cmip.bcc.cma.cn": 10,
#     "cmip.dess.tsinghua.edu.cn": 10,
#     "cmip.fio.org.cn": 10,
#     "dist.nmlab.snu.ac.kr": 10,
#     "esg-cccr.tropmet.res.in": 10,
#     "esg-dn1.nsc.liu.se": 10,
#     "esg-dn2.nsc.liu.se": 10,
#     "esg.camscma.cn": 10,
#     "esg.lasg.ac.cn": 10,
#     "esg1.umr-cnrm.fr": 10,
#     "esgf-cnr.hpc.cineca.it": 10,
#     "esgf-data2.diasjp.net": 10,
#     "esgf-data3.ceda.ac.uk": 10,
#     "esgf-data3.diasjp.net": 10,
#     "esgf-nimscmip6.apcc21.org": 10,
#     "esgf-node2.cmcc.it": 10,
#     "esgf.bsc.es": 10,
#     "esgf.dwd.de": 10,
#     "esgf.ichec.ie": 10,
#     "esgf.nci.org.au": 10,
#     "esgf.rcec.sinica.edu.tw": 10,
#     "esgf3.dkrz.de": 10,
#     "noresg.nird.sigma2.no": 10,
#     "polaris.pknu.ac.kr": 10,
#     "vesg.ipsl.upmc.fr": 10,
# }


import requests
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from requests.exceptions import ConnectionError

import pandas as pd
from typing import List, Dict, Tuple
from pangeo_forge_recipes.patterns import pattern_from_file_sequence
from pangeo_forge_recipes.recipes import XarrayZarrRecipe

def recipe_from_urls(urls, kwargs):
    # parse kwargs for different steps of the recipe
    pattern_kwargs = kwargs.get('pattern_kwargs', {})
    recipe_kwargs = kwargs.get('recipe_kwargs', {})
    
    pattern = pattern_from_file_sequence(urls, "time", **pattern_kwargs)
    recipe = XarrayZarrRecipe(
        pattern, xarray_concat_kwargs={"join": "exact"}, **recipe_kwargs
    )
    return recipe

def iid_request(session:requests.Session, netcdf3_session:requests.Session, iid:str, node:str) -> Tuple[List, Dict[str, Dict[str, str]]]:
    """Takes an instance id and returns urls and kwargs dict based on ESGF API request to the given `node`"""
    params = {
        "type": "File",
        "retracted": "false",
        "replica": "false",
        "format": "application/solr+json",
        # "fields": ["url", "size", "retracted", "table_id", "title","instance_id"], # TODO: why does this not work? Ill revisit when I am tuning performance, for now get all
        "latest": "true",
        "distrib": "true",
        # "limit": 500, # TODO: Should this be less?
    }
    facets = facets_from_iid(iid)
    params.update(facets)
    # searching with version does not work. So what I will do here is delete the version, and later check if the version of the iid is equal to the files found, otherwise error out. TODO
    del params['version']
    
    print('API request')
    resp = session.get(node, params=params)
    
    if not resp.status_code == 200:
        raise RuntimeError(f'Request failed with {resp.status_code} for {iid}')
        ## TODO: Implement more sophisticated error handling
    
    # TODO: Check header?
    resp_data = resp.json()['response']['docs']
    
    if len(resp_data) == 0:
        raise ValueError(f'No Files were found for {iid}')
        
        # Extract info
    raw_urls,sizes,retracted, table_ids, titles = zip(*[(rd['url'],rd['size'], rd['retracted'], rd['table_id'],rd['title']) for rd in resp_data])
        
    # Check consistency with iid input
    _check_response_facets_consistency(facets, resp_data)
    # this takes care of checking that all table_ids are the same, so I can do this
    table_id = table_ids[0][0]
    
    #pick http url
    urls = [_parse_url_type(url[0]) for url in raw_urls]
    
    # Check retractions (this seems a bit redundant, but what the heck
    if not all(r is False for r in retracted):
        print(retracted)
        raise ValueError(f"Query for {iid} contains retracted files")
    
    
    # extract date range from filename
    # TODO: Is there a more robust way to do this?
    # otherwise maybe use `id` (harder to parse)
    dates = [t.replace(".nc", "").split("_")[-1].split("-") for t in titles]


    # infer number of timesteps using pandas
    def format_date(str_date):
        return "-".join([str_date[0:4], str_date[4:]])

    # TODO: For non-monthly data, we have to make the freq input smarter
    timesteps = [
        len(pd.date_range(format_date(a[0]), format_date(a[1]), freq="1MS"))
        for a in dates
    ]
    print(f'Dates for each file: {dates}')
    print(f"Size per file in MB: {[f/1e6 for f in sizes]}")
    print(f"Inferred timesteps per file: {timesteps}")
    element_sizes = [size / n_t for size, n_t in zip(sizes, timesteps)]
    
    ### Determine kwargs
    print("Generate kwargs")
    # MAX_SUBSET_SIZE=1e9 # This is an option if the revised subsetting still runs into errors.
    MAX_SUBSET_SIZE=500e6
    DESIRED_CHUNKSIZE=200e6
    # TODO: We need a completely new logic branch which checks if the total size (sum(filesizes)) is smaller than a desired chunk
    target_chunks = {
        "time": choose_chunksize(
            allowed_divisors[table_id],
            DESIRED_CHUNKSIZE,
            element_sizes,
            timesteps,
            include_last=False,
        )
    }
    
    
    # dont even try subsetting if none of the files is too large
    if max(sizes)<=MAX_SUBSET_SIZE:
        subset_input = 0
    else:
        ## Determine subset_input parameters given the following constraints
        # - Needs to keep the subset size below MAX_SUBSET_SIZE
        # - (Not currently implemented) Resulting subsets should be evenly dividable by target_chunks (except for the last file, that can be odd). This might ultimately not be required once we figure out the locking issues. I cannot fulfill this right now with the dataset structure where often the first and last files have different number of timesteps than the 'middle' ones. 
        
        smallest_divisor = int(max(sizes)//MAX_SUBSET_SIZE+1)# need to subset at least with this to stay under required subset size
        subset_input = smallest_divisor
    
    recipe_kwargs = {"target_chunks": target_chunks}
    if subset_input > 1:
        recipe_kwargs["subset_inputs"] = {"time": subset_input}
    
    print(
        f"Will result in max chunksize of {max(element_sizes)*target_chunks['time']/1e6}MB"
    )
    
    # Check for netcdf version TODO: This is quite slow, not sure why...
    # print(urls)
    print(f"Check for netcdf 3 files")
    pattern_kwargs = {}
    if is_netcdf3(netcdf3_session, urls[-1]):
        pattern_kwargs['file_type']="netcdf3"
    
    print(f"Sort Urls by time")
    # sort urls in decending time order (to be able to pass them directly to the pangeo-forge recipe)
    end_dates = [a[-1] for a in dates]
    urls = [url for _, url in sorted(zip(end_dates, urls))]
    
    kwargs = {'recipe_kwargs':recipe_kwargs, 'pattern_kwargs':pattern_kwargs}
    print(f"Dynamically determined kwargs: {kwargs}")
    return urls, kwargs

def _parse_url_type(url:str) -> str:
    """Checks that url is of a desired type (currently only http) and removes appended text"""
    # From naomis code, in case we need to support OPENDAP
            #         resp = resp["docs"]
            #         offset += len(resp)
            #         # print(offset,numFound,len(resp))
            #         for d in resp:
            #             dataset_id = d["dataset_id"]
            #             dataset_size = d["size"]
            #             for f in d["url"]:
            #                 sp = f.split("|")
            #                 if sp[-1] == files_type:
            #                     url = sp[0]
            #                     if sp[-1] == "OPENDAP":
            #                         url = url.replace(".html", "")
            #                     dataset_url = url
            #             all_frames += [[dataset_id, dataset_url, dataset_size]]
    
    split_url = url.split('|')
    if not split_url[-1] == 'HTTPServer':
        raise ValueError('This recipe currently only supports HTTP links')
    else:
        return split_url[0]


def facets_from_iid(iid:str) -> Dict[str,str]:
    """Translates iid string to facet dict according to CMIP6 naming scheme"""
    iid_name_template = "mip_era.activity_id.institution_id.source_id.experiment_id.variant_label.table_id.variable_id.grid_label.version"
    facets = {}
    for name, value in zip(iid_name_template.split("."), iid.split(".")):
        facets[name] = value
    return facets

def choose_chunksize(
    chunksize_candidates: List[int],
    max_size: float,
    element_size_lst: List[float],
    timesteps_lst: List[int],
    include_last: bool = True,
) -> int:
    """Determines the ideal chunksize based on a list of preferred `divisors` and
    informations about the input files
    given the following constraints:
    - The resulting chunks are smaller than `max_size`
    - The determined chunksize will divide each file into even chunks
      (if `include_last` is false, the last file is allowed to have uneven chunks,
      but cannot be larger than the number of timesteps in the last file)

    Parameters
    ----------
    candidate_chunks : List[int]
        A list of chunksizes to consider.
    max_size : float
        Maximum size (in bytes) of the resulting chunksize
    element_size_lst : List[float]
        List of sizes (in bytes) of a single element along the chunking dimension (often time)
        for each of the input elements (files).
    timesteps_lst : List[int]
        List of timesteps for input elements
    include_last : bool, optional
        Option to include or exclude the last element from above lists, by default True.
        If number of elements of lists above is 1, this is always True

    Returns
    -------
    int
        Choosen chunksize
    """
    #     # TODO: infer clean divisions of the divisor (e.g. [1, 2, 3, 4, 6] for 12) automatically here
    #     candidate_chunks = divisors[:-1]+list(range(divisors[-1], max(timesteps_lst), divisors[-1]))

    if (
        not include_last and len(timesteps_lst) > 1
    ):  # we cannot exclude the last one if there is only one element.
        chunksize_filtered = [
            cs
            for cs in chunksize_candidates
            if all(
                nt % cs == 0 for nt in timesteps_lst[:-1]
            )  # do I need and timesteps_lst[-1] > cs
        ]
    else:
        chunksize_filtered = [
            cs
            for cs in chunksize_candidates
            if all(nt % cs == 0 for nt in timesteps_lst)
        ]
    output_chunksizes = [
        max([cs for cs in chunksize_filtered if cs * element_size <= max_size])
        for element_size in element_size_lst
    ]
    # what do we do if somehow this ends up being different? Take the min/max?
    if not all(oc == output_chunksizes[0] for oc in output_chunksizes):
        raise ValueError("Determined chunksizes are not all equal.")
    else:
        return output_chunksizes[0]

def is_netcdf3(session: requests.Session, url:str) -> bool:
    """Simple check to determine the netcdf file version behind a url.
    Requires the server to support range requests"""
    headers = {"Range": "bytes=0-2"}
    resp = session.get(url, headers=headers)
    return 'CDF' in str(resp.content)

def _check_response_facets_consistency(facets:Dict[str, str], file_resp:Dict[str, str]):
    # Check that all responses indeed have the same attributes
    # (error out on e.g. mixed versions for now)
    # TODO: We might allow mixed versions later, but need to be careful with that!
    check_facets = [
        "mip_era",
        "activity_id",
        "institution_id",
        "source_id",
        "experiment_id",
        "variant_label",
        "table_id",
        "variable_id",
        "grid_label",
        "version",
    ]

    def _check_single_element_list(lst):
        # double check that the facet returns are just a single element
        [out] = lst  # errors on a list with more than one element
        return out

    for fac in check_facets:
        file_facets = [_check_single_element_list(f[fac]) for f in file_resp]
        if not all(ff == file_facets[0] for ff in file_facets):
            raise ValueError(
                f"Found non-matching values for {fac} in search query response. Got {file_facets}"
            )  


## global variables
#TODO: Do we need more search nodes?
node_dict = {
    "llnl": "https://esgf-node.llnl.gov/esg-search/search",
    "dkrz": "https://esgf-data.dkrz.de/esg-search/search",
    "ipsl": "https://esgf-node.ipsl.upmc.fr/esg-search/search",
    "ceda": "https://esgf-index1.ceda.ac.uk/esg-search/search",
}

# For certain table_ids it is preferrable to have time chunks that are a multiple of e.g. 1 year for monthly data.
monthly_divisors = sorted(
    [1, 3, 6, 12, 12 * 3] + list(range(12 * 5, 12 * 200, 12 * 5)) + [684, 1026, 2052] 
    # the last list accomodates some special cases for `DAMIP` files (which are often only one file, but with a very odd number of years (e.g.  171 years for hist-aer 🤷). 
    #TODO: I might not want to allow this in the ocean and ice fields. Lets see
)

allowed_divisors = {
    "Omon": monthly_divisors,
    "SImon": monthly_divisors,
    "Amon": monthly_divisors,
}  # Add table_ids and allowed divisors as needed


## Recipe Generation
iids = [
    'CMIP6.DAMIP.BCC.BCC-CSM2-MR.hist-aer.r1i1p1f1.Amon.pr.gn.v20190507',
    'CMIP6.DAMIP.BCC.BCC-CSM2-MR.hist-aer.r2i1p1f1.Amon.pr.gn.v20190507',
    'CMIP6.DAMIP.BCC.BCC-CSM2-MR.hist-aer.r3i1p1f1.Amon.pr.gn.v20190508',
    'CMIP6.DAMIP.CAS.FGOALS-g3.hist-aer.r1i1p1f1.Amon.pr.gn.v20200411',
    'CMIP6.DAMIP.CAS.FGOALS-g3.hist-aer.r2i1p1f1.Amon.pr.gn.v20200411',
    'CMIP6.DAMIP.CAS.FGOALS-g3.hist-aer.r3i1p1f1.Amon.pr.gn.v20200411',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r10i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r10i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r11i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r11i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r12i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r12i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r13i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r13i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r14i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r14i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r15i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r15i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r1i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r1i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r2i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r2i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r3i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r3i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r4i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r4i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r5i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r5i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r6i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r6i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r7i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r7i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r8i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r8i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r9i1p1f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CCCma.CanESM5.hist-aer.r9i1p2f1.Amon.pr.gn.v20190429',
    'CMIP6.DAMIP.CNRM-CERFACS.CNRM-CM6-1.hist-aer.r10i1p1f2.Amon.pr.gr.v20190308',
    'CMIP6.DAMIP.CNRM-CERFACS.CNRM-CM6-1.hist-aer.r1i1p1f2.Amon.pr.gr.v20190308',
    'CMIP6.DAMIP.CNRM-CERFACS.CNRM-CM6-1.hist-aer.r2i1p1f2.Amon.pr.gr.v20190308',
    'CMIP6.DAMIP.CNRM-CERFACS.CNRM-CM6-1.hist-aer.r3i1p1f2.Amon.pr.gr.v20190308',
    'CMIP6.DAMIP.CNRM-CERFACS.CNRM-CM6-1.hist-aer.r4i1p1f2.Amon.pr.gr.v20190308',
    'CMIP6.DAMIP.CNRM-CERFACS.CNRM-CM6-1.hist-aer.r5i1p1f2.Amon.pr.gr.v20190308',
    'CMIP6.DAMIP.CNRM-CERFACS.CNRM-CM6-1.hist-aer.r6i1p1f2.Amon.pr.gr.v20190308',
    'CMIP6.DAMIP.CNRM-CERFACS.CNRM-CM6-1.hist-aer.r7i1p1f2.Amon.pr.gr.v20190308',
    'CMIP6.DAMIP.CNRM-CERFACS.CNRM-CM6-1.hist-aer.r8i1p1f2.Amon.pr.gr.v20190308',
    'CMIP6.DAMIP.CNRM-CERFACS.CNRM-CM6-1.hist-aer.r9i1p1f2.Amon.pr.gr.v20190308',
    'CMIP6.DAMIP.CSIRO-ARCCSS.ACCESS-CM2.hist-aer.r1i1p1f1.Amon.pr.gn.v20201120',
    'CMIP6.DAMIP.CSIRO-ARCCSS.ACCESS-CM2.hist-aer.r2i1p1f1.Amon.pr.gn.v20201120',
    'CMIP6.DAMIP.CSIRO-ARCCSS.ACCESS-CM2.hist-aer.r3i1p1f1.Amon.pr.gn.v20201120',
    'CMIP6.DAMIP.CSIRO.ACCESS-ESM1-5.hist-aer.r1i1p1f1.Amon.pr.gn.v20200615',
    'CMIP6.DAMIP.CSIRO.ACCESS-ESM1-5.hist-aer.r2i1p1f1.Amon.pr.gn.v20200615',
    'CMIP6.DAMIP.CSIRO.ACCESS-ESM1-5.hist-aer.r3i1p1f1.Amon.pr.gn.v20200615',
    'CMIP6.DAMIP.IPSL.IPSL-CM6A-LR.hist-aer.r10i1p1f1.Amon.pr.gr.v20180914',
    'CMIP6.DAMIP.IPSL.IPSL-CM6A-LR.hist-aer.r1i1p1f1.Amon.pr.gr.v20180914',
    'CMIP6.DAMIP.IPSL.IPSL-CM6A-LR.hist-aer.r2i1p1f1.Amon.pr.gr.v20180914',
    'CMIP6.DAMIP.IPSL.IPSL-CM6A-LR.hist-aer.r3i1p1f1.Amon.pr.gr.v20180914',
    'CMIP6.DAMIP.IPSL.IPSL-CM6A-LR.hist-aer.r4i1p1f1.Amon.pr.gr.v20180914',
    'CMIP6.DAMIP.IPSL.IPSL-CM6A-LR.hist-aer.r5i1p1f1.Amon.pr.gr.v20180914',
    'CMIP6.DAMIP.IPSL.IPSL-CM6A-LR.hist-aer.r6i1p1f1.Amon.pr.gr.v20180914',
    'CMIP6.DAMIP.IPSL.IPSL-CM6A-LR.hist-aer.r7i1p1f1.Amon.pr.gr.v20180914',
    'CMIP6.DAMIP.IPSL.IPSL-CM6A-LR.hist-aer.r8i1p1f1.Amon.pr.gr.v20180914',
    'CMIP6.DAMIP.IPSL.IPSL-CM6A-LR.hist-aer.r9i1p1f1.Amon.pr.gr.v20180914',
    'CMIP6.DAMIP.MIROC.MIROC6.hist-aer.r10i1p1f1.Amon.pr.gn.v20201228',
    'CMIP6.DAMIP.MIROC.MIROC6.hist-aer.r1i1p1f1.Amon.pr.gn.v20190705',
    'CMIP6.DAMIP.MIROC.MIROC6.hist-aer.r2i1p1f1.Amon.pr.gn.v20190705',
    'CMIP6.DAMIP.MIROC.MIROC6.hist-aer.r3i1p1f1.Amon.pr.gn.v20190705',
    'CMIP6.DAMIP.MIROC.MIROC6.hist-aer.r4i1p1f1.Amon.pr.gn.v20201228',
    'CMIP6.DAMIP.MIROC.MIROC6.hist-aer.r5i1p1f1.Amon.pr.gn.v20201228',
    'CMIP6.DAMIP.MIROC.MIROC6.hist-aer.r6i1p1f1.Amon.pr.gn.v20201228',
    'CMIP6.DAMIP.MIROC.MIROC6.hist-aer.r7i1p1f1.Amon.pr.gn.v20201228',
    'CMIP6.DAMIP.MIROC.MIROC6.hist-aer.r8i1p1f1.Amon.pr.gn.v20201228',
    'CMIP6.DAMIP.MIROC.MIROC6.hist-aer.r9i1p1f1.Amon.pr.gn.v20201228',
    'CMIP6.DAMIP.MOHC.HadGEM3-GC31-LL.hist-aer.r1i1p1f3.Amon.pr.gn.v20190814',
    'CMIP6.DAMIP.MOHC.HadGEM3-GC31-LL.hist-aer.r2i1p1f3.Amon.pr.gn.v20190815',
    'CMIP6.DAMIP.MOHC.HadGEM3-GC31-LL.hist-aer.r3i1p1f3.Amon.pr.gn.v20190814',
    'CMIP6.DAMIP.MOHC.HadGEM3-GC31-LL.hist-aer.r4i1p1f3.Amon.pr.gn.v20190814',
    'CMIP6.DAMIP.MOHC.HadGEM3-GC31-LL.hist-aer.r5i1p1f3.Amon.pr.gn.v20211123',
    'CMIP6.DAMIP.MRI.MRI-ESM2-0.hist-aer.r1i1p1f1.Amon.pr.gn.v20190320',
    'CMIP6.DAMIP.MRI.MRI-ESM2-0.hist-aer.r2i1p1f1.Amon.pr.gn.v20200327',
    'CMIP6.DAMIP.MRI.MRI-ESM2-0.hist-aer.r3i1p1f1.Amon.pr.gn.v20190320',
    'CMIP6.DAMIP.MRI.MRI-ESM2-0.hist-aer.r4i1p1f1.Amon.pr.gn.v20200327',
    'CMIP6.DAMIP.MRI.MRI-ESM2-0.hist-aer.r5i1p1f1.Amon.pr.gn.v20190320',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r1i1p1f1.Amon.pr.gn.v20180821',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r1i1p1f2.Amon.pr.gn.v20191226',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r1i1p3f1.Amon.pr.gn.v20191226',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r2i1p1f1.Amon.pr.gn.v20180821',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r2i1p1f2.Amon.pr.gn.v20191226',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r2i1p3f1.Amon.pr.gn.v20191226',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r3i1p1f1.Amon.pr.gn.v20180822',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r3i1p1f2.Amon.pr.gn.v20191226',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r3i1p3f1.Amon.pr.gn.v20191226',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r4i1p1f1.Amon.pr.gn.v20180823',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r4i1p1f2.Amon.pr.gn.v20191226',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r4i1p3f1.Amon.pr.gn.v20191226',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r5i1p1f1.Amon.pr.gn.v20180823',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r5i1p1f2.Amon.pr.gn.v20191226',
    'CMIP6.DAMIP.NASA-GISS.GISS-E2-1-G.hist-aer.r5i1p3f1.Amon.pr.gn.v20191226',
    'CMIP6.DAMIP.NCAR.CESM2.hist-aer.r1i1p1f1.Amon.pr.gn.v20200206',
    'CMIP6.DAMIP.NCAR.CESM2.hist-aer.r3i1p1f1.Amon.pr.gn.v20200305',
    'CMIP6.DAMIP.NCC.NorESM2-LM.hist-aer.r1i1p1f1.Amon.pr.gn.v20190920',
    'CMIP6.DAMIP.NCC.NorESM2-LM.hist-aer.r2i1p1f1.Amon.pr.gn.v20190920',  
    'CMIP6.DAMIP.NCC.NorESM2-LM.hist-aer.r3i1p1f1.Amon.pr.gn.v20190920',
    'CMIP6.DAMIP.NOAA-GFDL.GFDL-ESM4.hist-aer.r1i1p1f1.Amon.pr.gr1.v20180701',
]

# set up a requests session to speed up the recipe generation
# Some tips: https://stackoverflow.com/questions/62599036/python-requests-is-slow-and-takes-very-long-to-complete-http-or-https-request
session = requests.Session()
netcdf3_session = requests.Session()

# TODO: The range requests are taking FOREVER. I need to speed those up.

recipes = {}
for ii, iid in enumerate(iids):
    success = False
    print(f"\n+++ Generate Recipe for {iid} ({ii+1}/{len(iids)}) +++")
    for node_name, node_url in node_dict.items():
        print(f"Node: {node_name}")
        try:
            urls, kwargs = iid_request(session, netcdf3_session, iid, node_url)
            if urls:
                recipes[iid] = recipe_from_urls(urls, kwargs)
                success = True
                # if successful break the node loop.
                break
            else:
                print(f"No urls provided for {iid} and Node:{node_name}")
        except ConnectionError as e: # I cant for the hell of it find out what the exception is...
            print(e)
            print(f'It seems that there was a connection issue for {iid} and Node:{node_name}')
    if not success:
        print(f'!!!!!!!!!!!!!! Recipe creation failed for {iid}!!!!!!!!!!!!!!!!!\nn')