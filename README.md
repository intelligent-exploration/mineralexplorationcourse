# GEOL3888 Economic Geology
This repository comprises the work required for Prac classes from Week 6, 7, 8, 9, and 10.
Modern day exploration requires a versatile skill-set, including being able to analyse, visualise and interpret data. This section of the course should give you hard-skills in the Python programming language. Plus you will learn about Machine Learning, Plate Tectonic Reconstructions, all while using specific Python packages (numpy, scipy, scikit-learn, matplotlib, cartopy, pygplates, etc), and a few other data and research tools and environments (e.g. Git, GPlates, databases, Docker containers).

The prac classes are compulsory. And the assessment will be a write-up of one of the pracs. More details given in the assessment section.

# Environment installation and setup

Download this repo which contains all the small datasets and the codes/notebooks we will be using: https://github.com/intelligent-exploration/mineralexplorationcourse/archive/master.zip

There are a few additional datasets that will be needed for some of the examples and exercises. Download them as instructed in the notes, but mostly they can all be found at https://www.dropbox.com/sh/mhvss5491nij1sb/AABC4VkGNkU4A1Ej5gusk3BFa?dl=0

## Set up the python environment via Docker
Pull it from [docker hub](https://cloud.docker.com/u/nbutter/repository/docker/nbutter/pyforgeo) or build the [Dockerfile](Docker_details/Dockerfile) and run!

or....

## Create the python environment locally
We need two different environments. One is particularly needed for pyGplates. We are using the conda package manager. Conda is installed on the computers in the lab. If you wish you use your own computer you can download it from https://repo.continuum.io/miniconda/ 

### python 3 environment (week 8 and 9)
Once Conda is installed, to create the environment to work in use:

```
conda create -n pyGEOL python=3.7 scipy=1.2 scikit-learn=0.20 matplotlib=3.0 pyshp=2.0 numpy=1.16 jupyter=1.0 cartopy=0.17 pandas=0.24 notebook=5.7.4 git
```

Then to activate and run the notebooks, invoke:
```
conda activate pyGEOL
jupyter notebook
```

### python 2 pyGplates environment (week 10)

***python 2 specifically is required for the week 10 notebook only because of the pygplates dependency.***

Note: Windows users before creating your conda environment you will need to set
```
set CONDA_FORCE_32BIT=1
```

Make the conda environement:
```
conda create -n py2GEOL python=2.7 scipy=1.2 scikit-learn=0.20 matplotlib=2.0 pyshp=1.2 numpy=1.15 jupyter=1.0 cartopy=0.17 pandas=0.24 notebook=5.7.4 git
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

### class

 * [Class exercise and examples ](Week6/exercise)


# Week 7

Python lists, 2d arrays and numpy arrays
 * https://www.machinelearningplus.com/python/numpy-tutorial-part1-array-python-examples/
 * https://www.tutorialspoint.com/numpy/index.htm
 * https://matplotlib.org/tutorials/introductory/pyplot.html
 * https://docs.scipy.org/doc/numpy/reference/generated/numpy.savetxt.html 
 * https://www.ict.social/python/basics/multidimensional-lists-in-python
 * https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html

 * [Class exercise and examples ](Week7/exercise)

# Week 8

## Intro to ML

* General Intro to ML: https://www.seas.upenn.edu/~cis519/fall2017/lectures/01_introduction.pdf https://scikit-learn.org/stable/tutorial/statistical_inference/supervised_learning.html

* Linear regression: https://www.datacamp.com/community/tutorials/essentials-linear-regression-python https://www.cs.toronto.edu/~frossard/post/linear_regression/  https://towardsdatascience.com/linear-regression-using-gradient-descent-97a6c8700931
https://www.youtube.com/watch?v=zPG4NjIkCjc https://www.youtube.com/watch?v=JvS2triCgOY Mean squared error https://www.youtube.com/watch?v=r-txC-dpI-E

* Gradient descent: https://en.wikipedia.org/wiki/Gradient_descent https://ml-cheatsheet.readthedocs.io/en/latest/gradient_descent.html

* Logistic regression: https://towardsdatascience.com/logistic-regression-detailed-overview-46c4da4303bc
* https://www.kaggle.com/emilyhorsman/basic-logistic-regression-with-numpy https://medium.com/@martinpella/logistic-regression-from-scratch-in-python-124c5636b8ac https://medium.com/technology-nineleaps/logistic-regression-gradient-descent-optimization-part-1-ed320325a67e https://www.youtube.com/watch?v=yIYKR4sgzI8 https://www.youtube.com/watch?v=yIYKR4sgzI8&t=58s https://www.cc.gatech.edu/~bboots3/CS4641-Fall2016/Lectures/Lecture6.pdf

* https://www.geeksforgeeks.org/learning-model-building-scikit-learn-python-machine-learning-library/

## Python for GIS

A Python introduction leading into visualising geo-spatial data.
* [An Introduction to Python for GIS](Week8/Intro_Python_Geo.ipynb)

Then a deeper lesson exploring some new Python libraries and advanced features of Python data manipulation on various types of data.
* [Python with shapefiles and pandas](Week8/PandasExamples.ipynb)


# Week 9


We will use Python to perform machine learning on a well-constructed dataset 
* [Machine Learning with Python for Geoscience](Week9/ML_Geo.ipynb)
* Random forest background: https://medium.com/@williamkoehrsen/random-forest-simple-explanation-377895a60d2d 
* Random forest code: https://www.datascience.com/resources/notebooks/random-forest-intro

# Week 10

This is an in-depth exploration linking some of the earlier workflows. These notebooks and workflows will teach you how to to create and manipulate data, use python from outside of Jupyter (i.e. normal python), and to use bespoke Python packages (pygplates and build-your-own). This task can be expanded to your full project assessment if chosen.

* [Porphyry Copper Deposits in the Andes](Week10/)

# Week 11

* Github tutorial: https://product.hubspot.com/blog/git-and-github-tutorial-for-beginners
* Video tutorial: https://www.youtube.com/watch?v=0fKg7e37bQE&t=523s
* Gitkraken Software tutorial: https://www.youtube.com/watch?v=ub9GfRziCtU

# Week 13

* https://www.waikato.ac.nz/library/study/guides/write-scientific-reports
* https://writingcenter.unc.edu/tips-and-tools/scientific-reports/
* http://www.lcc.uma.es/~eat/pdf/sw.pdf
* [Evaluation criteria for presentation](Week13/)


# Assessment
1) Pass/Fail. The 3 notebooks from week 8 and 9  contain a specific assessment at the end. Complete the task and show your tutor. 

2) 10%. Choose **one (1)** of the datasets and key figures from the workbooks of week 8 and 9 to write up as a full paper-style report (e.g. Intro_Python_Geo.ipynb, PandasExamples.ipynb, OR ML_Geo.ipynb). This should include a short abstract, a methods sections, and an interpretation of the results and any insights you can make about the data. You can do this all within the Jupyter notebook environment - approx 4 pages max INCLUDING figures and NOT counting code (but please include your code, especially any comments you make to exisiting code).  **Due: the start of the Week 10 lecture, i.e. Thursday 9th May 2019 10am**

3) 25% (15% Technical Report, 10% Software/code/methods). Group Assignment. All members will recieve the same mark. The technical report should feature an abstract, introduction, a brief background of case study, and results and discussion.  Follow standard research paper presentation format with literature review and citations. See [Butterworth et al. 2016](https://doi.org/10.1002/2016TC004289) for inspiration. Project details to be given. Possible project ideas will be running the Week 10 workflow on different data sets. Discuss your project with your tutor. 
**Due: start of Week 13 Prac, i.e Thursday 30th May 2019 1pm**

3) 15%. Group Presentation. Presentation of project. We want to know about the project you did! Max 20 minute presentation. With questions to follow. All members are encouraged to present (your presentation mark will be based on the strongest speaker, so this is great practice for presenting your work). **Held in Week 13 prac time, i.e Thursday 30th May 2019 1pm-4pm.**



