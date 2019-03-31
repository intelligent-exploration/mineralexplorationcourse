This section of the Prac will be updated with new content soon... here you can see the beta version.

# AndesML


This repo is an example of how machine learning can be used to find Porphyry Copper deposits in the Andes. Modified from the workflows of Butterworth, N., D. Steinberg, R. D. Müller, S. Williams, A. S. Merdith, and S. Hardy (2016), Tectonic environments of South American porphyry copper magmatism through time revealed by spatiotemporal data mining, Tectonics, 35, 2847–2862, doi:10.1002/2016TC004289


## Exercises
We will be running a workflow to build up a dataset for visualising. 

## Setup
Step 1: Download the required files from here.
In the Week 9 Folder should be 

### Datafiles
convergence.py
CopperDeposits
Muller_Bertrand_coregistered.pkl
Muller_Bertrand_coregistered_random.pkl
Muller_Bertrand_coregistered_sampleMuller0.pkl
Muller_convergence
Muller_copper_prob.ipynb
Muller_etal_2016_AREPS_Agegrids_v1.11
Muller_export
Muller_gplates
README.md
Utils_coreg.py
Utils_coreg.pyc


### Python environemnt setup
Create the required conda environemnt:

conda install python=2.7 
conda install -y scipy=1.2 scikit-learn=0.20 matplotlib=2.0 pyshp=1.2 numpy=1.15 jupyter=1.0 cartopy=0.17 pandas=0.24 notebook=5.7.4


Install pygplates beta-revision-12
https://sourceforge.net/projects/gplates/files/pygplates/beta-revision-12/


You can also use the Docker container if that is easier!

Now let's start the actual prac!
sudo docker run  -it --rm -v`pwd`:/workspace pyg /bin/bash -c "source activate py2GEOL && convergence.py Muller_gplates/Global_EarthByte_230-0Ma_GK07_AREPS.rot Muller_export/topology_subduction_boundaries_sL_1.00Ma.shp Muller_export/topology_subduction_boundaries_sR_1.00Ma.shp Muller_export/topology_platepolygons_1.00Ma.shp 1"



Run the analysis with the Seton plate model.
Run with a different plate model
Run with correct age grids
Run with EarthChem data
Run with Bertrand data in SE Asia

Include paper with notes.

## Instructions 
The main notebook, Muller_copper_prob.ipynb, can be used with the pre-made datasets.

***Muller_Bertrand_coregistered.pkl*** contains the full set of copper-deposits and their associated tectono-magmatic properites.

***Muller_Bertrand_coregistered_random.pkl*** contains a psudeo-random set of non-deposits with known tectono-magmatic properties that can be used for training.

***Muller_Bertrand_coregistered_sampleMuller0.pkl*** contains a set of points following the present day South American subduction zone and the correspoiding tectono-magmatic properites of those points. 

The full set of data developed here (and age-grids) can be downloaded directly from https://www.dropbox.com/s/gm7ydn376jarkkf/Muller.zip?dl=0

Or you can make your own files by following the instructions in this readme and other notebooks/scripts in this repo. 

A python environemt to get this to work can be made using

```
set CONDA_FORCE_32BIT=1
conda create -n pyforml_32 python=2.7
conda activate pyforml_32
conda install jupyter matplotlib scipy pyshp numpy scikit-learn=0.17
```
Then make sure you install pygplates, on Windows I unzipped the 32 bit version then include this line before importing:
```
import sys
sys.path.append("C:\pygplates_rev12_python27_win32")
import pygplates
```

## Instructions if you want to do it your self.

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

### Run the convergence script on the topologies for all time steps, using something like this.
```
for age in {0..230}; do echo ${age};  python convergence.py Muller_gplates/Global_EarthByte_230-0Ma_GK07_AREPS.rot Muller_export/topology_subduction_boundaries_sL_${age}.00Ma.shp Muller_export/topology_subduction_boundaries_sR_${age}.00Ma.shp Muller_export/topology_platepolygons_${age}.00Ma.shp ${age}; done
```

This should output a collection of files labelled as **"subStats_XXX.csv"** for each time interval XXX. These are the calculated plate kinematics data.

##What are the values in each of the columns representive of?
You might have to read the convergence.py script to see what it does. Check the header information or maybe somewhere in the script that prints out the file to try and reverse engineer it.

The *.pkl files, as produced in Butterworth et al 2016 require an additional step to create. This is the bit that 'co-registers' the copper deposits with the kinematic data from convergece.py and the agegrids. You have to run the "coregLoop" function within Utils_coreg.py 

There are a few manual changes that may need to be done as the Agegrid data, and directories are hard coded. But in principle it is easy , just use the lines in python: 

```
#Read in a shapefile
andeanPoints="XYBer14_t2_ANDES.shp"
[recs,shapes,fields,Nshp]=readTopologyPlatepolygonFile(andeanPoints)

#Run the coregistration
coregData=coregLoop(recs,shapes,fields,Nshp,3)

#Then you can save the data as a pickle file or do whatever you like with it!
```
