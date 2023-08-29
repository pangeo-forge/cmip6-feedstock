# 🚨🚨🚨Deprecation Notice🚨🚨🚨
This feedstock is no longer maintained. All further efforts of the [Pangeo / ESGF Cloud Data Working Group](https://pangeo-data.github.io/pangeo-cmip6-cloud/) to ingest CMIP6 data into Analysis-Ready Cloud-Optimized zarr stores are consolidated [here](https://github.com/leap-stc/cmip6-leap-feedstock).


# cmip6-feedstock

**This feedstock will not actually contain recipes anymore, and only be used for all requests to add CMIP data to the cloud that is not currently available in the [Pangeo CMIP6 Cloud holdings](https://pangeo-data.github.io/pangeo-cmip6-cloud/).**

## How to request new datasets
Ideally a request for new datasets includes a list of [instance_ids]() that can be added directly to one of the feedstocks below. [pangeo-forge-esgf](https://github.com/jbusecke/pangeo-forge-esgf) provides a convenient way of getting a list of available datasets directly from ESGF, using a wildcard notation ([example](https://github.com/jbusecke/pangeo-forge-esgf#parsing-a-list-of-instance-ids-using-wildcards)).


## Split feedstocks
[CMIP6-PMIP](https://github.com/pangeo-forge/CMIP6-PMIP-feedstock)
[CMIP6-DAMIP](https://github.com/pangeo-forge/CMIP6-DAMIP-feedstock)


Please refer to https://pangeo-forge.org/dashboard/feedstock/7 to follow this feedstock's journey within Pangeo Forge Cloud.

Pangeo Forge Feedstocks are collaborative.

If you see a way that the recipes within this feedstock can be improved, please open an Issue or Pull Request on this repo.
