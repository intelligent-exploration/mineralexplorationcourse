#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
20190315 nbutter
Collection of tools and functions.
Use 

from utils import *

and call functions directly.

'''
print("importing modules...")

import shapefile
print("Imported shapefile.")

import matplotlib.colors as mcolors
print("Imported matplotlib.")

print("defining functions...")


def readTopologyPlatepolygonFile(filename):
    '''
    Reads shapefiles and returns the all the data fields
    '''
    shapeRead = shapefile.Reader(filename)

    recs    = shapeRead.records()
    shapes  = shapeRead.shapes()
    fields  = shapeRead.fields
    Nshp    = len(shapes)
    
    return(recs,shapes,fields,Nshp)


def colormap_age():
	'''
	Make your own colormap!
	A map made from printing out the 

	for time in xrange(0,20,1):
	    [(time+0.1)/20.0,1-(time+0.1)/20.0,(time+0.1)/20.0]

	or from the functions that read a gmt cpt file

	Inspired by http://stackoverflow.com/a/11659600
	'''

	levs = range(20)
	assert len(levs) % 2 == 0, 'N levels must be even.'

	age_cmap = mcolors.LinearSegmentedColormap.from_list(name='green_purple', 
	                                                 colors =[(0.0, 0.0, 0.0) ,
	(0.82352941176470584, 0.0, 0.0) ,
	(0.90196078431372551, 0.15686274509803921, 0.0) ,
	(0.96078431372549022, 0.23529411764705882, 0.0) ,
	(1.0, 0.3843137254901961, 0.0) ,
	(1.0, 0.58039215686274515, 0.0) ,
	(1.0, 0.77647058823529413, 0.0) ,
	(0.97254901960784312, 0.94901960784313721, 0.0) ,
	(0.88627450980392153, 0.94901960784313721, 0.0) ,
	(0.63921568627450975, 1.0, 0.0) ,
	(0.63921568627450975, 1.0, 0.0) ,
	(0.34509803921568627, 1.0, 0.0) ,
	(0.34509803921568627, 1.0, 0.0) ,
	(0.0, 0.94117647058823528, 0.0) ,
	(0.0, 0.94117647058823528, 0.0) ,
	(0.0, 1.0, 0.82745098039215681) ,
	(0.0, 0.88235294117647056, 1.0) ,
	(0.0, 0.58823529411764708, 1.0) ,
	(0.0, 0.59999999999999998, 1.0) ,
	(0.0, 0.40000000000000002, 1.0) ,
	(0.0, 0.20000000000000001, 1.0) ,
	(0.0, 0.0, 1.0) ,
	(0.20000000000000001, 0.0, 1.0) ,
	(0.40000000000000002, 0.0, 1.0) ,
	(0.59999999999999998, 0.0, 1.0) ,
	(0.80000000000000004, 0.0, 1.0) ,
	(0.80000000000000004, 0.0, 1.0) ],N=len(levs)-1,)

	return(age_cmap)

print("Finished utils.py load!")