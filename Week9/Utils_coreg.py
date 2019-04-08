#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
20140902 nbutter
Collection of tools and functions.

Updated May 2018

'''

#import sys
#sys.path.append("../PYGPLATES/pygplates_rev18_python27_win32")
import pygplates
print("Imported pygplates.")

import shapefile
print("Imported shapefile.")

import numpy
print("Imported numpy.")

import scipy
import scipy.io
import scipy.spatial
import scipy.stats
print("Imported scipy.")
 

### ---------------- Utility functions ---------------- ###
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

def printsom():
        print('Hello')

def readCSV(filename):
    '''
    Reads point data
    '''
    #lonMid, latMid, convRate, distance, \
    #orthAbs, orthOP, orthTrench,subObliquity,subPolarity,distEdge = \
    f=numpy.loadtxt(filename, delimiter=',')
    
    
    #return(lonMid, latMid, convRate, distance, \
    #orthAbs, orthOP, orthTrench,subObliquity,subPolarity,distEdge)
    return(f)
    
##
def gridRead(filename):
    '''
    #Reads a netcdf grid file expecting the grid variables named 'lat', 'lon', 
        'z'.
    
    input:
        string filename of a grid
    returns: 
        three arrays describing each variable in file
    usage:
        [Lat,Lon,Data]=gridRead('~/mydir/myfile.grd')
        
    '''
    #import scipy.io #For reading netcdf grid files

    #Use netcdf module to load in the data
    #A netcdf module exists also, but is pretty much the same but deals with 
    #newer version also.
    data = scipy.io.netcdf.netcdf_file(filename,'r')
    
    
    #Get out each grid axis name
    #var=[]
    #for k, v in data.variables.iteritems():
    #    var.append(k)
    #print "Variables should be [Lat, Long, Data]"
    #print "Variables in "+ filename +" are ", var, "and look like:"
    ##Store the data in varibales
    #varX=data.variables[str(var[0])][:]
    #varY=data.variables[str(var[1])][:]
    #varZ=numpy.array(data.variables[str(var[2])][:])
    #print "Lat: ",varX[0:4],"to",varX[-5:-1]
    #print "Lon: ",varY[0:4],"to",varY[-5:-1]
    #print "Data: ",varZ[0]
    
    try:
        varX=data.variables['lon'][:]
        #Make sure all are -180 to 180
        varX=convert360to180(varX)
        varY=data.variables['lat'][:]
        varZ=numpy.array(data.variables['z'][:])

    except:
        print("Warning:grid data does not contain variables 'lat', 'lon' or 'z'")
        print("variables names are:")
        for k, v in data.variables.iteritems():
            print(k)
        print("Returning dummy variables for", filename)
        varX=numpy.arange(0,360,1)
        varY=numpy.arange(-90,90,1)
        varZ=numpy.zeros([180,360])
        
    return(varX,varY,varZ)
    

def coregPoint(point,data,region):
    '''
    Finds the nearest neighbour to a point from a bunch of other points
    point - array([longitude,latitude])
    data - array
    region - integer, same units as data
    '''
    tree = scipy.spatial.cKDTree(data)
    dists, indexes = tree.query(point,k=1,distance_upper_bound=region) 

    if indexes==len(data):
        return 'inf'
    else:
        return (indexes)

#
def points_in_circle(circle, arr):
    '''
    A generator to return all points whose indices are within given circle.
    http://stackoverflow.com/a/2774284
    Warning: If a point is near the the edges of the raster it will not loop 
    around to the other side of the raster!
    '''
    i0,j0,r = circle

    for i in xrange(intceil(i0-r),intceil(i0+r)):
        ri = numpy.sqrt(r**2-(i-i0)**2)
        for j in xrange(intceil(j0-ri),intceil(j0+ri)):
            if (i >= 0 and i < len(arr[:,0])) and (j>=0 and j < len(arr[0,:])):
                yield arr[i][j]

#            
def intceil(x):
    return int(numpy.ceil(x))                                            

#
def coregRaster(point,data,region):
    '''
    Finds the mean value of a raster, around a point with a specified radius.
    point - array([longitude,latitude])
    data - array
    region - integer, same units as data
    '''
    i0=point[1]
    j0=point[0]
    r=region #In units of degrees
    pts_iterator = points_in_circle((i0,j0,region), data)
    pts = numpy.array(list(pts_iterator))

    #return(scipy.stats.nanmean(pts)) #deprecated from scipy 0.15
    return(numpy.nanmean(pts))
###
def regressData(data):
    '''
    TODO....
    http://stackoverflow.com/a/11479279
    
    #Remove nans from data
    f=g[~np.isnan(g).any(1)]
    #Combine the rows we want
    regData=numpy.concatenate((f[:,2:4],f[:,6:9]),axis=1)
    #Then run this function...
    #Coeffiencts are this:
    clf.coef_
    '''
    from sklearn import linear_model
    clf = linear_model.LinearRegression()
    clf.fit(regData,f[:,4])
           
###
def makeTemplate(data):
    """
    A point distribution model to use as a template for procrustes analysis, or an average of the data
    """
    tlist=range(0,data.shape[1])
    combinedData=numpy.mean(data,0)
    a=numpy.dstack((tlist,combinedData))
    
    return(a.reshape((data.shape[1],2)))


##
#A function to filter data based on criteria in a column
def cleanCondition(cleanArray,dataArray):
    #Pick the points we want using the conditonal array we just built
    andesClean=dataArray[cleanArray,:,:]
    #Now convert the tuple we made into an array
    andesClean=numpy.asarray(andesClean)
    #And, clean it up
    andesClean=numpy.reshape(andesClean,[andesClean.shape[1],andesClean.shape[2],andesClean.shape[3]])
    return andesClean


def convert180to360(longData):
    '''
    Converts Longitude data from '-180 to 180' to '0 to 360'
    Keep in mind the data will not be reordered in any way.

    input:
        Array of longitudes
    returns: 
        Array of longitudes with values replaced
    usage:
        [LonFix]=convert180to360(Lons)

    '''
    arrayFix=[]
    for index,item in enumerate(longData):
        value=360+item
        if (item < 0):
             arrayFix.append(value)
        else:
             arrayFix.append(item)
    arrayFix=numpy.array(arrayFix)
    return(arrayFix)

def convert360to180(longData):
    '''
    Converts Longitude data from '-180 to 180' to '0 to 360'
    Keep in mind the data will not be reordered in any way.

    input:
        Array of longitudes
    returns: 
        Array of longitudes with values replaced
    usage:
        [LonFix]=convert360to180(Lons)

    '''
    arrayFix=[]
    for index,item in enumerate(longData):
        value=item-360
        if (item > 180):
             arrayFix.append(value)
        else:
             arrayFix.append(item)
    arrayFix=numpy.array(arrayFix)
    return(arrayFix)


def colormap_age():
	'''
	Make your own colormap!
	A map made from printing out the 

	for time in xrange(0,20,1):
	    [(time+0.1)/20.0,1-(time+0.1)/20.0,(time+0.1)/20.0]

	or from the functions that read a gmt cpt file

	Inspired by http://stackoverflow.com/a/11659600
	'''
	import matplotlib.colors as mcolors

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


