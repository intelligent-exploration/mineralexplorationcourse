# GEOL3888 Economic Geology
This repository comprises the work required for Prac classes from Week 6, 7, 8, and 9.
Modern day exploration requires a versatile skill-set, including being able to analyse, visualise and interpret data. This section of the course should give you hard-skills in the Python programming language. Plus you will learn about Machine Learning, Plate Tectonic Reconstructions, all while using specific Python packages (numpy, scipy, scikit-learn, matplotlib, cartopy, pygplates, etc), and a few other data and research tools and environments (e.g. Git, GPlates, databases, Docker containers).

The prac classes are compulsory. And the assessment will be a write-up of one of the pracs. More details given in the assessment section.

# Environment installation and setup

Download this repo: ht which contains all the small datasets and the codes/notebooks we will be using: https://github.com/intelligent-exploration/mineralexplorationcourse/archive/master.zip

Next, download the bigger datasets that are needed for some of the examples and exercises:
***UPLOAD data to dropbox/host?***

## Set up the python environment via Docker
Pull it from [docker hub](https://cloud.docker.com/u/nbutter/repository/docker/nbutter/pyforgeo) or build the [Dockerfile](Docker_details/Dockerfile) and run!

or....

## Create the python environment locally
We need two different environments. One is particularly needed for pyGplates. We are using the conda package manager. Conda is installed on the computers in the lab. If you wish you use your own computer you can download it from https://repo.continuum.io/miniconda/ 

### python 3 environment (week 7 and 8)
Once Conda is installed, to create the environment to work in use:

```
conda create -n pyGEOL scipy=1.2 scikit-learn=0.20 matplotlib=3.0 pyshp=2.0 numpy=1.16 jupyter=1.0 cartopy=0.17 pandas=0.24 notebook=5.7.4
```

Then to activate and run the notebooks, invoke:
```
conda activate pyGEOL
jupyter notebook
```

### python 2 pyGplates environment (week 9)

***python 2 specifically is required for the week 9 notebook only because of the pygplates dependency.***

Note: Windows users before creating your conda environment you will need to set
```
set CONDA_FORCE_32BIT=1
```

Make the conda environement:
```
conda create -n py2GEOL python=2.7 scipy=1.2 scikit-learn=0.20 matplotlib=2.0 pyshp=1.2 numpy=1.15 jupyter=1.0 cartopy=0.17 pandas=0.24 notebook=5.7.4
```

Download and unzip pygplates:

Windows
https://sourceforge.net/projects/gplates/files/pygplates/beta-revision-12/pygplates_rev12_python27_win32.zip

Mac
https://sourceforge.net/projects/gplates/files/pygplates/beta-revision-12/pygplates_rev12_python27_MacOS64.zip

Add pygplates to your Python Path in whatever environement you are using 

From within Python via e.g.:
```
import sys
sys.path.append("C:\pygplates_rev12_python27_win32")
import pygplates
```

Or from outside of Python
```
export PYTHONPATH=${PYTHONPATH}:/Users/Downloads/pygplates_rev12_python27_MacOS64/
```




# Week 6
This will be a hands-on programming intro. You will learn the basics of coding and Python. If you have experience using command line or some programming language (e.g. Matlab, C++) then you should find this straight-forward. If you have never programmed before, great! This will be the beginning of a great relationship with getting data to do what you want!

### Getting started videos

Mac
https://www.youtube.com/watch?v=ftJoIN_OADc
https://www.youtube.com/watch?v=SPZamL2Xbsc

Linux
https://www.youtube.com/watch?v=3xp-ixFbDuE

Windows
https://www.youtube.com/watch?v=dX2-V2BocqQ

### Online programming tutorials and examples
 
https://snakify.org/en/lessons/print_input_numbers/
https://www.geeksforgeeks.org/python-programming-examples/

### online python compiler 
https://repl.it/repls


# Week 7
A Python introduction leading into visualising geo-spatial data.
* [An Introduction to Python for GIS](Week7/Intro_Python_Geo.ipynb)

Then a deeper lesson exploring some new Python libraries and advanced features of Python data manipulation on various types of data.
* [Python with shapefiles and pandas](Week7/PandasExamples.ipynb)

# Week 8
We will use Python to perform machine learning on a well-constructed dataset 
* [Machine Learning with Python for Geoscience](Week8/ML_Geo.ipynb)

# Week 9
This is an in-depth exploration linking some of the earlier workflows. This notebook will teach you how to to create and manipulate data, use python from outside of Jupyter (i.e. normal python), and to use bespoke Python packages (pygplates and build-your-own). This task can be expanded to your full project assessment if chosen.

# Assessment
1) Pass/Fail. The 3 notebooks from week 7 and 8 each contain a specific assessment at the end. Complete the task and show your tutor. 

2) 10%. Choose one of the datasets and key figures from the workbooks of week 7 and 8, to write up as a full paper-style report. This should include a short abstract, a methods sections, and an interpretation of the results. You can do this all within the Jupyter notebook environment. This is due at the start of Week 10 prac, i.e. Thursday May 9, 2019 1pm.
