############################# just copied from esfg.py (TODO: being able to import this, would be nice and clean) ############################################
"""ESGF API Search Results to Pandas Dataframes
"""
import requests
import numpy
import pandas as pd
import warnings
from typing import List, Dict

# dummy comment
# copied from Naomis code https://github.com/pangeo-data/pangeo-cmip6-cloud/blob/master/myconfig.py
target_keys = [
    "activity_id",
    "institution_id",
    "source_id",
    "experiment_id",
    "member_id",
    "table_id",
    "variable_id",
    "grid_label",
]

target_format = "%(" + ")s/%(".join(target_keys) + ")s"


node_pref = {
    "esgf-data1.llnl.gov": 0,
    "esgf-data2.llnl.gov": 0,
    "aims3.llnl.gov": 0,
    "esgdata.gfdl.noaa.gov": 10,
    "esgf-data.ucar.edu": 10,
    "dpesgf03.nccs.nasa.gov": 5,
    "crd-esgf-drc.ec.gc.ca": 6,
    "cmip.bcc.cma.cn": 10,
    "cmip.dess.tsinghua.edu.cn": 10,
    "cmip.fio.org.cn": 10,
    "dist.nmlab.snu.ac.kr": 10,
    "esg-cccr.tropmet.res.in": 10,
    "esg-dn1.nsc.liu.se": 10,
    "esg-dn2.nsc.liu.se": 10,
    "esg.camscma.cn": 10,
    "esg.lasg.ac.cn": 10,
    "esg1.umr-cnrm.fr": 10,
    "esgf-cnr.hpc.cineca.it": 10,
    "esgf-data2.diasjp.net": 10,
    "esgf-data3.ceda.ac.uk": 10,
    "esgf-data3.diasjp.net": 10,
    "esgf-nimscmip6.apcc21.org": 10,
    "esgf-node2.cmcc.it": 10,
    "esgf.bsc.es": 10,
    "esgf.dwd.de": 10,
    "esgf.ichec.ie": 10,
    "esgf.nci.org.au": 10,
    "esgf.rcec.sinica.edu.tw": 10,
    "esgf3.dkrz.de": 10,
    "noresg.nird.sigma2.no": 10,
    "polaris.pknu.ac.kr": 10,
    "vesg.ipsl.upmc.fr": 10,
}


# Author: Unknown
# I got the original version from a word document published by ESGF
# https://docs.google.com/document/d/1pxz1Kd3JHfFp8vR2JCVBfApbsHmbUQQstifhGNdc6U0/edit?usp=sharing
# API AT:
# https://github.com/ESGF/esgf.github.io/wiki/ESGF_Search_REST_API#results-pagination


def esgf_search(
    search,
    server="https://esgf-node.llnl.gov/esg-search/search",
    files_type="HTTPServer",
    local_node=True,
    project="CMIP6",
    page_size=500,
    verbose=False,
    format="application%2Fsolr%2Bjson",
    toFilter=True,
):

    client = requests.session()
    payload = search
    payload["project"] = project
    payload["type"] = "File"
    if local_node:
        payload["distrib"] = "false"

    payload["format"] = format
    payload["limit"] = 500

    numFound = 10000
    all_frames = []
    offset = 0
    while offset < numFound:
        payload["offset"] = offset
        url_keys = []
        for k in payload:
            url_keys += ["{}={}".format(k, payload[k])]

        url = "{}/?{}".format(server, "&".join(url_keys))
        print(url)
        r = client.get(url)
        r.raise_for_status()
        resp = r.json()["response"]
        numFound = int(resp["numFound"])

        resp = resp["docs"]
        offset += len(resp)
        # print(offset,numFound,len(resp))
        for d in resp:
            dataset_id = d["dataset_id"]
            dataset_size = d["size"]
            for f in d["url"]:
                sp = f.split("|")
                if sp[-1] == files_type:
                    url = sp[0]
                    if sp[-1] == "OPENDAP":
                        url = url.replace(".html", "")
                    dataset_url = url
            all_frames += [[dataset_id, dataset_url, dataset_size]]

    ddict = {}
    item = 0
    for item, alist in enumerate(all_frames):
        dataset_id = alist[0]
        dataset_url = alist[1]
        dataset_size = alist[2]
        vlist = dataset_id.split("|")[0].split(".")[-9:]
        vlist += [dataset_url.split("/")[-1]]
        vlist += [dataset_size]
        vlist += [dataset_url]
        vlist += [dataset_id.split("|")[-1]]
        ddict[item] = vlist
        item += 1

    dz = pd.DataFrame.from_dict(ddict, orient="index")
    if len(dz) == 0:
        print("empty search response")
        return dz

    dz = dz.rename(
        columns={
            0: "activity_id",
            1: "institution_id",
            2: "source_id",
            3: "experiment_id",
            4: "member_id",
            5: "table_id",
            6: "variable_id",
            7: "grid_label",
            8: "version_id",
            9: "ncfile",
            10: "file_size",
            11: "url",
            12: "data_node",
        }
    )

    dz["ds_dir"] = dz.apply(lambda row: target_format % row, axis=1)
    dz["node_order"] = [node_pref[s] for s in dz.data_node]
    dz["start"] = [s.split("_")[-1].split("-")[0] for s in dz.ncfile]
    dz["stop"] = [s.split("_")[-1].split("-")[-1].split(".")[0] for s in dz.ncfile]

    if toFilter:
        # remove all 999 nodes
        dz = dz[dz.node_order != 999]

        # keep only best node
        dz = dz.sort_values(by=["node_order"])
        dz = dz.drop_duplicates(subset=["ds_dir", "ncfile", "version_id"], keep="first")

        # keep only most recent version from best node
        dz = dz.sort_values(by=["version_id"])
        dz = dz.drop_duplicates(subset=["ds_dir", "ncfile"], keep="last")

    return dz


####################################################
from pangeo_forge_recipes.patterns import pattern_from_file_sequence
from pangeo_forge_recipes.recipes import XarrayZarrRecipe

# from esgf import (
#     esgf_search,
# )  # We probably want to strip this out later, left as is for now.

node_dict = {
    "llnl": "https://esgf-node.llnl.gov/esg-search/search",
    "ipsl": "https://esgf-node.ipsl.upmc.fr/esg-search/search",
    "ceda": "https://esgf-index1.ceda.ac.uk/esg-search/search",
    "dkrz": "https://esgf-data.dkrz.de/esg-search/search",
}


def urls_from_instance_id(instance_id):
    # get facets from instance_id
    facet_labels = (
        "mip_era",
        "activity_id",
        "institution_id",
        "source_id",
        "experiment_id",
        "member_id",
        "table_id",
        "variable_id",
        "grid_label",
        "version",
    )

    facet_vals = instance_id.split(".")
    if len(facet_vals) != 10:
        raise ValueError(
            "Please specify a query of the form {"
            + ("}.{".join(facet_labels).upper())
            + "}"
        )

    facets = dict(zip(facet_labels, facet_vals))

    if facets["mip_era"] != "CMIP6":
        raise ValueError("Only CMIP6 mip_era supported")

    # version doesn't work here
    keep_facets = (
        "activity_id",
        "institution_id",
        "source_id",
        "experiment_id",
        "member_id",
        "table_id",
        "variable_id",
        "grid_label",
    )
    search_facets = {f: facets[f] for f in keep_facets}

    search_node = "llnl"
    ESGF_site = node_dict[
        search_node
    ]  # TODO: We might have to be more clever here and search through different nodes. For later.

    df = esgf_search(search_facets, server=ESGF_site)  # this modifies the dict inside?

    # get list of urls
    urls = df["url"].tolist()

    # sort urls in decending time order (to be able to pass them directly to the pangeo-forge recipe)
    end_dates = [url.split("-")[-1].replace(".nc", "") for url in urls]
    urls = [url for _, url in sorted(zip(end_dates, urls))]

    # version is still not working
    # if facets["version"].startswith("v"):
    #    facets["version"] = facets["version"][1:]

    # TODO Check that there are no gaps or duplicates.

    return urls


## Misc logic
def facets_from_iid(iid):
    iid_name_template = "mip_era.activity_id.institution_id.source_id.experiment_id.variant_label.table_id.variable_id.grid_label.version"
    facets = {}
    for name, value in zip(iid_name_template.split("."), iid.split(".")):
        facets[name] = value
    return facets


## Logic to dynamically generate input kwargs
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


def dynamic_kwarg_generation(
    iid: str,
) -> Dict[str, Dict[str, int]]:
    """Dynamically generates keyword arguments `target_chunks` and `subsset_input` for
    recipe generation based on information available via the ESGF API

    Parameters
    ----------
    iid : str
        ESGF instance_id

    Returns
    -------
    Dict[str,Dict[str, int]]
        Dictionary containing keyword arguments that can be passed to `XarrayZarrRecipe`

    """
    print(f"====== Generate kwargs for {iid} =======")
    # TODO, I query the API in multiple places. Need to refactor something robust (which might try different urls?)

    url = "https://esgf-node.llnl.gov/esg-search/search"
    # url = "https://esgf-data.dkrz.de/esg-search/search"

    # TODO: the 'distrib' parameter does not work as expected for all datasets.
    # Need to investigate that on the ESGF side.
    # Could just iterate through nodes for now.
    
    # MAX_SUBSET_SIZE=1e9 # This is an option if the revised subsetting still runs into errors.
    MAX_SUBSET_SIZE=500e6
    DESIRED_CHUNKSIZE=200e6

    params = {
        "type": "File",
        "retracted": "false",
        "replica": "false",
        "format": "application/solr+json",
        # "fields": "size",
        "latest": "true",
        # "distrib": "true",
        "limit": 500,
    }

    facets = facets_from_iid(iid)
    params.update(facets)

    del params[
        "version"
    ]  # TODO: Why do we have to delete this? Need to understand that better
    resp = requests.get(url=url, params=params)

    file_resp = resp.json()["response"]["docs"]

    if not len(file_resp) > 0:
        raise ValueError("ESGF API query did not return any files.")

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

    # now make sure that the table_id is a key in `preferred_time_divisions` otherwise error
    table_id = file_resp[0]["table_id"][
        0
    ]  # Confirmed before that this list is only 1 element
    if table_id not in allowed_divisors.keys():
        raise ValueError(
            f"Didnt find `table_id` value {table_id} in the `allowed_divisors` dict."
        )

    filesizes = [f["size"] for f in file_resp]

    # extract date range from filename
    # TODO: Is there a more robust way to do this?
    # otherwise maybe use `id` (harder to parse)
    dates = [a["title"].replace(".nc", "").split("_")[-1].split("-") for a in file_resp]


    # infer number of timesteps using pandas
    def format_date(str_date):
        return "-".join([str_date[0:4], str_date[4:]])

    # TODO: For non-monthly data, we have to make the freq input smarter
    timesteps = [
        len(pd.date_range(format_date(a[0]), format_date(a[1]), freq="1MS"))
        for a in dates
    ]
    print(f'Dates for each file: {dates}')
    print(f"Size per file in MB: {[f/1e6 for f in filesizes]}")
    print(f"Inferred timesteps per file: {timesteps}")
    element_sizes = [size / n_t for size, n_t in zip(filesizes, timesteps)]
    
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
    if max(filesizes)<=MAX_SUBSET_SIZE:
        subset_input = 0
    else:
        ## Determine subset_input parameters given the following constraints
        # - Needs to keep the subset size below MAX_SUBSET_SIZE
        # - (Not currently implemented) Resulting subsets should be evenly dividable by target_chunks (except for the last file, that can be odd). This might ultimately not be required once we figure out the locking issues. I cannot fulfill this right now with the dataset structure where often the first and last files have different number of timesteps than the 'middle' ones. 
        
        smallest_divisor = int(max(filesizes)//MAX_SUBSET_SIZE+1)# need to subset at least with this to stay under required subset size
        subset_input = smallest_divisor
        

    dynamic_kwargs = {"target_chunks": target_chunks}
    if subset_input > 1:
        dynamic_kwargs["subset_inputs"] = {"time": subset_input}
    print(f"Dynamically determined kwargs: {dynamic_kwargs} for {iid}")
    print(
        f"Will result in max chunksize of {max(element_sizes)*target_chunks['time']/1e6}MB"
    )
    return dynamic_kwargs


## global variables

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
    # These dont produce a valid url. TODO: we need to find a good way to handle these (as they might be flaky)
    #'CMIP6.DAMIP.NCC.NorESM2-LM.hist-aer.r2i1p1f1.Amon.pr.gn.v20190920',  
    #'CMIP6.DAMIP.NCC.NorESM2-LM.hist-aer.r3i1p1f1.Amon.pr.gn.v20190920',
    'CMIP6.DAMIP.NOAA-GFDL.GFDL-ESM4.hist-aer.r1i1p1f1.Amon.pr.gr1.v20180701',
]
inputs = {iid: dynamic_kwarg_generation(iid) for iid in iids}


def recipe_from_urls(urls, instance_kwargs):
    pattern = pattern_from_file_sequence(urls, "time")

    recipe = XarrayZarrRecipe(
        pattern, xarray_concat_kwargs={"join": "exact"}, **instance_kwargs
    )
    return recipe


recipes = {
    iid: recipe_from_urls(urls_from_instance_id(iid), kwargs)
    for iid, kwargs in inputs.items()
}
