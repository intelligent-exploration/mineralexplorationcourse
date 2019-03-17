Australian Groundwater Explorer River Region Data Download
--------------------------------------------------------------------------

BACKGROUND
------------
The Bureau of Meteorology collates groundwater data from state agencies in accordance with the Water Regulations and publishes it as nationally-consistent datasets including:
- National Groundwater Information System (NGIS)
- water levels and salinity measurements

NGIS is a national, spatially-enabled groundwater database. It contains the following data:
- bores
- lithology logs
- construction logs
- hydrostratigraphy (borehole) logs
- hydrogeologic unit table
- 3D construction lines
- 3D hydrostratigraphy (borehole) lines
- georasters (rasters defining the top/base of hydrogeologic units)
- geovolumes (3D geometry of hydrogeologic units)
- management zones (groundwater management areas)
- aquifers (aquifer boundaries)

The Bureau also harvests hydrochemistry data from Geoscience Australia (GA) web services including:
- bore information e.g. location, purpose
- analyses e.g. major and minor ions, environmental tracers
Note: GA uses project-based bore IDs, which may be different to the bore IDs used by State/Territory agencies and the NGIS.

Further information about the NGIS, water level, salinity and hydrochemistry data is available from the Explorer Metadata web page:http://www.bom.gov.au/water/groundwater/explorer/metadata.shtml

DATASET
-----------
The current dataset is an extract of NGIS, water levels, salinity and hydrochemisty data for a particular river region, which has been processed for download from the Australian Groundwater Explorer:
http://www.bom.gov.au/water/groundwater/explorer/map.shtml

Note river region boundaries are taken from the Australian Hydrological Geospatial Fabric (Geofabric):
http://www.bom.gov.au/water/geofabric/index.shtml


Zip file: gw_<format>_<river region>.zip

This file contains:
- NGIS version 1.5 extract
- Water level measurements:level_<river region>.csv
- Salinity measurements: salinity_<river region>.csv
- GA hydrochemistry bore locations and attributes: ga_bore_<river region>.csv
- GA hydrochemistry measurements: ga_hydrochem_<river region>.csv
- Product release notes: gw_region_README.txt (this file)

The Australian Groundwater Explorer provides data in the following formats: 
- ESRI file geodatabase + csvs
- shapefiles + csvs
- KML files (bores only)

Only a subset of NGIS is available to download from the Australian Groundwater Explorer and the available datasets depend on the format as listed above. 

For geodatabase format, the following NGIS data is available:
- bores
- lithology logs
- construction logs
- hydrostratigraphy (borehole) logs
- hydrogeologic unit table
- 3D construction lines
- 3D hydrostratigraphy (borehole) lines
- management zones
- aquifers

For shapefile format, the following NGIS data is available:
- bores
- lithology logs
- construction logs
- hydrostratigraphy (bore) logs
- hydrogeologic unit table
- management zones
- aquifers

Only bore data is available in KML format.

NGIS georasters and geovolumes are not included in river region or state download package. These datasets are part of the national NGIS database, which is available by email request to the Bureau's groundwater mailbox: groundwater@bom.gov.au


Water level and salinity measurements are provided in csvs with geodatabase and shapefile format downloads. The currency of water levels and salinity data depends on the river region. For both measurements, in many regions, the most recent data is available in July 2018. It should be noted that these data is not available for some regions.


Groundwater level information are presented using following different terms or variables:

- Depth to water (DTW)—measured from the ground surface to the groundwater level. Positive values are below the ground surface. Negative values are above the ground surface and indicate artesian conditions. 
- Standing water level (SWL)—measured from the reference point on the bore (e.g. the top of casing) to the groundwater level. Positive values are below the reference point and negative values are above the     reference point.
- Reduced standing water level (RSWL)—the groundwater level elevation relative to the Australian Height Datum (AHD). This is the standard height datum used in Australia. It sets mean sea level as zero elevation. Positive values are above mean sea level and negative values are below.


Salinity data is measured in one of two ways:

- Total Dissolved Solids (TDS) in units of milligrams per litre (mg/L), which measures the weight of dissolved content per litre of water.
- Electrical Conductivity (EC) in units of micro-Siemens or milli-Siemens per centimetre (µS/cm or mS/cm). Electrical Conductivity is a reliable estimator of the salinity of water.


Hydrochemistry data provided by the Geoscience Australia is also available in csvs with geodatabase and shapefile format downloads. This data is provided in two csv files per region where one file gives bore location and attributes while the other provdes hydrochemistry data.

The downloads from the Groundwater Explorer are available by state and State/Territory. At this stage the download packages only contain data from the lead water agency in each state/territory and GA. This means that they do not include data from NT Power and Water, 
National Collaborative Research Infrastructure Strategy (NCRIS), WA Water Corporation and WA Department of Primary Industries and Regional Development. Data for these agencies is available by email request to the Bureau's groundwater mailbox: 
groundwater@bom.gov.au

The complete, national NGIS database is available by request to the Bureau's groundwater mailbox: 
groundwater@bom.gov.au


LINKING NGIS TO LEVEL AND SALINITY DATA
-----------------------------------------
Water levels and salinity data can be linked to the NGIS bore data to get the location and attributes of bores. 
The common identifier called hydroid in the water levels and salinity datasets and HydroID in the NGIS Bore datataset can be used for this purpose.


HYDROCHEMISTRY DATA
---------------------
The location and attributes of the bores with GA hydrochemistry data can be found in the ga_bore_<river region>.csv.
This can be linked to the GA hydrochemistry data using a common identifier called ga_bore_id.


INSTALLATION
--------------
Unzip the downloaded file, making sure to maintain the internal directory structure.


COPYRIGHT
--------
The data in the NGIS (bores, logs, aquifers and groundwater management areas) are subject to copyright. Refer to the NGIS Copyright page for more information: http://www.bom.gov.au/water/groundwater/ngis/copyright.shtml 

The water level, salinity and hydrochemistry data is also subject to copyright. Refer to the Explorer Copyright page for more information: http://www.bom.gov.au/water/groundwater/explorer/copyright.shtml


DOCUMENTATION/FURTHER INFORMATION
-----------------------------------
Further information is available from the Explorer About web page:
http://www.bom.gov.au/water/groundwater/explorer/index.shtml

NGIS documentation is available from the NGIS web page:
http://www.bom.gov.au/water/groundwater/ngis/documentation.shtml

Documentation includes:
- Data model diagram
- Data dictionary 
- Data product specification


CHANGE LOG FOR NGIS
---------------------

	CHANGE AT VERSION 1.5 (current version)
	-----------------
	
	- Updated NGIS bore and bore log data for all states and territories
	- Added data availability attributes to NGIS_Bore
	- Added hydrostratigraphy logs for New South Wales Murray Basin 
	- Added georasters and geovolumes for the Great Artesian Basin and St Vincent Basin
	- Updated georasters for the Daly Basin
	- Updated groundwater management areas 
	- Corrected issues with NT bore elevations
	- Updated metadata
	
	CHANGES AT VERSION 1.4.1 
	-----------------	
	- Updated Victorian Construction Log data as bore screens were missing
        - Added data from new data providers (NCRIS, WA Department of Primary Industries and Regional Development)

	CHANGES AT VERSION 1.4
	-----------------
	- Added aquifer boundary dataset (not available in river region download)
	- Updated NGIS bore and bore log data for all states and territories
	- Updated groundwater management areas including updating the data model for NGIS_ManagementZones dataset (not available in river region download)
	- Updated GA hydrochemistry data
	- Updated metadata
	
	CHANGES AT VERSION 1.3 
	-----------------
	- Added new bore and bore log data for all states and territories
	- Added groundwater management areas including updating the data model for NGIS_ManagementZones dataset (not available in river region download)
	- Updated metadata 
	- Added unique feature identifiers (HydroIDs) for all layers
	- Updated National Aquifer Framework (NAF) attributes to NAF version 1.1.2
	- Addressed errors identified during QA/QC process (unreasonable values, terms not standardised, etc)

	CHANGES AT VERSION 1.2
	-----------------
	- Updated metadata 
	- Added unique feature identifiers (HydroIDs) for all layers
	- Updated National Aquifer Framework (NAF) attributes to NAF version 1.1
	- Addressed errors identified during QA/QC process (unreasonable values, terms not standardised, etc)
	- Added georasters and geovolumes for hydrogeologic complexes in the Murray Basin
	- Included the version 2.0 of the Victorian Aquifer Framework rasters in NGIS_Georaster dataset
	- Added new bore logs for the Murray Basin
	- Added production bore data for the ACT
	- Minor update to the data model for NGIS_ManagementZones dataaset (no data included in current release)
	- Updated the State bore identifiers for Victoria (HydroCode and StateBoreID) to be consistent
	with the new Victorian water management system
	- Changed the State bore identifiers for Western Australia (HydroCode and StateBoreID) to be the AWRC 
	Reference Number 
		
	CHANGES AT VERSION 1.1
	-----------------
	- Metadata updated
	- Added Georaster catalog dataset with rasters for all Victoria 
	- ACT and WA Department of Water data added for the first time, all bore and bore log data.
	- Standardised state identifiers for all layers
	- Upgraded to NAF version 1.0 from previous draft NAF
	- Addressed errors identified during QA/QC process (unreasonable values, terms not standardised)

KNOWN ISSUES
---------------
- The data model for NGIS_Bore only accomodates information about one hydrogeologic unit aquifer. However, some bores screen multiple aquifers.
- HydroIDs may not be persistent or meaningful. They serve only as unique identifier in the NGIS database. They are initially populated by State/Territory lead water agencies. The Bureau adds a sate/Territory offset to ensure that they are nationally unique.
- The identifiers used in NGIS (HydroCode, StateBoreID, HydroID) may not be the same as the identifiers used in time series databases e.g. for Tasmania, different identifiers are used in the 
bore/bore log and time series databases. This can make it difficult to relate bore data to time series data.
- The feature of interest in NGIS is bore pipes. For multi-pipe bores, this has resulted in some repetition of data e.g. multiple lithology logs for the same hole.
- NGIS data was sourced from State/Territory lead water agencies in late 2017 - early 2018. Consequently NGIS does not hold recently added or updated data. This data must be sourced directly from lead water agencies.
- In WA, there may be overlap in the data supplied by Department of Water and Environmental Regulation and Water Corporation.
- Known issues with shapefile download format:
    - National NGIS is a file geodatabase and some shortcomings are introduced when the data is converted into shapefile format:
      -	Coded Value Domains (CVDs) do not work and the numeric codes appear in the fields with CVDs. Please refer to the NGIS Data Dictionary to see the corresponding values: ftp://ftp.bom.gov.au/anon/home/water/ngis/NgisDataDictionary.pdf 
    - Field names are truncated to a maximum length of 10 characters. Where truncation would make the field name the same as another field name, an underscore and sequential number is added to the end of the field name.



USER FEEDBACK
-------------
Your feedback is vital to assist us in maintaining a business case to support and enhance the Explorer and NGIS. Also by letting us know how you are using, or want to use the Explorer and NGIS, you will enable us to better serve user needs. If you have any ideas, issues or questions related to the Explore, NGIS or their use, please contact us through our email groundwater@bom.gov.au


CONTACT INFORMATION
-------------------
Grundwater Unit
Bureau of Meteorology
GPO Box 1289
Melbourne VIC 3001
Web: http://www.bom.gov.au/water/groundwater/
Email: groundwater@bom.gov.au

