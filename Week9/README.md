# Machine learning to classify ore deposits from tectonmagmatic properties

This repo is an example of how machine learning can be used to find Porphyry Copper deposits in the Andes. Modified from the workflows of Butterworth, N., D. Steinberg, R. D. Müller, S. Williams, A. S. Merdith, and S. Hardy (2016), Tectonic environments of South American porphyry copper magmatism through time revealed by spatiotemporal data mining, Tectonics, 35, 2847–2862, doi:10.1002/2016TC004289


## Setup
Step 1: Download this repo and additional data files from [here](https://www.dropbox.com/sh/vcgddw8tkh8lp51/AABBLNcKvIvZwAVTEoKRBuGYa?dl=0).

Unzip the download into the Week 9 folder. The Week 9 Folder should contain the following folders/files:

### Data folders
**Muller_etal_2016_AREPS_Agegrids_v1.11** - Age grids used as a coregistration dataset

**Muller_export** - Output data from convergence.py used as a coregistration dataset

**Muller_gplates** - Plate model topologies and rotation files.

**CopperDeposits** - A shapefile of copper deposits to be used as a coregistration dataset

### Python codes/scripts
**convergence.py** - python script that calculates convergence rates for a particular plate model. 

**coregloop.py** - python script to take a shapefile and "coregister" the points with some other datasets

**Muller_copper_prob.ipynb** - a python noebook that steps through the machine learning process.

**Utils_coreg.py** - Some additional tools and functions called in the main scripts.

All are these are in this Repository in the /Week9 folder.


## Python environemnt 

Create and activate the required conda environemnt:

```
conda create -n py2GEOL -y python=2.7 scipy=1.2 scikit-learn=0.20 matplotlib=2.0 pyshp=1.2 numpy=1.15 jupyter=1.0 cartopy=0.17 pandas=0.24 notebook=5.7.4

conda activate py2GEOL
```

Install pygplates beta-revision-12
https://sourceforge.net/projects/gplates/files/pygplates/beta-revision-12/



Alternatively you can also use the [Docker container](https://hub.docker.com/r/nbutter/pyforgeo) if that is easier!
```
docker pull nbutter/pyforgeo
```

If using docker most commands below can be completed by prefixing the comands with this docker line:
```
sudo docker run  -it --rm -v`pwd`:/workspace pyg /bin/bash -c "source activate py2GEOL && COMMAND-GOES-HERE"
```



## Step 1, get plate kinematic data

Run the ***convergence.py*** script from a bash terminal shell. This script calculates convergence rates and other kinemtics about plate model along the subduction boundaries. To run this script you need to specify a rotation file, left- and right-handed subduction zone topology files, a closed plate polygon file, and a time. We want to calculate kinematic data for the entire history of the plate model, so we can do this in a small bash script:

```
cd Week9/

for age in {0..230}; do echo ${age};  python convergence.py Muller_gplates/Global_EarthByte_230-0Ma_GK07_AREPS.rot Muller_export/topology_subduction_boundaries_sL_${age}.00Ma.shp Muller_export/topology_subduction_boundaries_sR_${age}.00Ma.shp Muller_export/topology_platepolygons_${age}.00Ma.shp ${age}; done
```

This will create a bunch of files in called ***subStats_XX.csv*** that contain the statistcs of the subduction zones (i.e. the kinemtatic data). See the comments in ***convergence.py*** to see exactly what is in output. These will be dumped in the top-level directory, so let's move them to a folder:
```
mv subStats*.csv Muller_convergence/
```

## Step 2, coregister the kinematic data with ore deposit data

The folder ***CopperDeposits*** contains a shapefile ***XYBer14_t2_ANDES.shp*** with a set up known copper deposits in the Andes. These are currently just points on a map with a formation age associated with them. Our ultimate goal is to identify what kind of tectonomagmatic environemnts are associated with the formation of these ore deposits, to begin to do this we first must link age-coded tectonomagmatic properties with our age-coded ore deposits. We can use the script ***coregLoop.py*** to achieve this!

Run the script coregLoop.py

Make sure all the paths to the various files are correct in the script. Then run with:
```
python coregLoop.py
```

This creates three datasets (in the python 'pickle' file format) that we can have enough data in them to perform some machine learning on:


***Muller_Bertrand_coregistered.pkl*** contains the full set of copper-deposits and their associated tectono-magmatic properites.

***Muller_Bertrand_coregistered_random.pkl*** contains a psudeo-random set of non-deposits with known tectono-magmatic properties that can be used for training.

***Muller_Bertrand_coregistered_sampleMuller0.pkl*** contains a set of points following the present day South American subduction zone and the correspoiding tectono-magmatic properites of those points. 


## Step 3, apply machine learning algorithms to find significant associations

Open up Muller_copper_prob.ipynb in a jupyter notebook and continue the intereactive analysis in there.



## Additional instructions to create different datasets.

This analysis can be done with any combination of rotation files plate polygons. To create the initial set use these instructions.

### Load the Muller et al. 2016 rotation and topology files into GPlates and Export subduction topologies, as outlined in "convergence.py":

* Load in plate polygons and rotation file to GPlates. Click on the export data button and select the following.
* 1. Export Resolved Topologies (CitcomS specific)
* 2. Shapefiles
* 3. ONLY select "Export all Plate Polygons to a single file"
* and "Export Plate Polygon segments to files based on feature type..."

### Copy out all the line segments of the data you have just exported. 
This hack is required because in some cases GPlates (the older versions did not do this) will split the topologies into "points" and "lines", we just want have to rename the "lines" files if this occurs for "convergence.py" to read. This assumes you have exported 230 Myr of reconstructions. Use this bash snippet or similar to find and rename the shape files.

```
for age in {0..230}; do echo ${age};cp Muller_export/topology_subduction_boundaries_sL_${age}.00Ma/topology_subduction_boundaries_sL_${age}.00Ma_polyline.shx Muller_export/topology_subduction_boundaries_sL_${age}.00Ma.shx;
cp Muller_export/topology_subduction_boundaries_sL_${age}.00Ma/topology_subduction_boundaries_sL_${age}.00Ma_polyline.shp Muller_export/topology_subduction_boundaries_sL_${age}.00Ma.shp;
cp Muller_export/topology_subduction_boundaries_sL_${age}.00Ma/topology_subduction_boundaries_sL_${age}.00Ma_polyline.prj Muller_export/topology_subduction_boundaries_sL_${age}.00Ma.prj;
cp Muller_export/topology_subduction_boundaries_sL_${age}.00Ma/topology_subduction_boundaries_sL_${age}.00Ma_polyline.dbf Muller_export/topology_subduction_boundaries_sL_${age}.00Ma.dbf; done
```
```
for age in {0..230}; do echo ${age};cp Muller_export/topology_subduction_boundaries_sR_${age}.00Ma/topology_subduction_boundaries_sR_${age}.00Ma_polyline.shx Muller_export/topology_subduction_boundaries_sR_${age}.00Ma.shx;
cp Muller_export/topology_subduction_boundaries_sR_${age}.00Ma/topology_subduction_boundaries_sR_${age}.00Ma_polyline.shp Muller_export/topology_subduction_boundaries_sR_${age}.00Ma.shp;
cp Muller_export/topology_subduction_boundaries_sR_${age}.00Ma/topology_subduction_boundaries_sR_${age}.00Ma_polyline.prj Muller_export/topology_subduction_boundaries_sR_${age}.00Ma.prj;
cp Muller_export/topology_subduction_boundaries_sR_${age}.00Ma/topology_subduction_boundaries_sR_${age}.00Ma_polyline.dbf Muller_export/topology_subduction_boundaries_sR_${age}.00Ma.dbf; done
```

You are then ready to run ***convergence.py*** as in the instructions above.

### Create an alternate ore deposit dataset.

The Ore Deposits require a their lat/lon location, along with an age, and to be used in this workflow need to be assigned a GPlates "Plate ID", the example set already has this. To assign a plate id to a new dataset you can use pygplates or do this within GPlates.

