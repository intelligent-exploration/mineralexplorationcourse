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

def coregLoop(recs,shapes,fields,Nshp,shapeType):
    '''
    coregLoop
    
    This script reconstructs a shapefile of points to their birth time and 
    coregisters each point with another set of points, or a raster file.
    
    INPUTS:
    Can best be used with the output of readTopologyPlatepolygonFile() 
    which reads a shapefile and returns the 'records', 'shapes', 'fields',
    and 'number of shapes'.
 
    Hardcoded filenames/variables must be changed below:
        Rotation file - input_rotation_filename
        Time dependent raster files - rasterfile
        Time dependent vector csv outputs from 'convergence.py' - f
    
    OUTPUTS:
    Coregistered array: List of lat/lons with properties.

    METHOD:
    Takes a set of points (from a shapefile)
    Rotates the points back to their birth position
    Determines the point's birth position geophysical properties (coregisters)
    '''
    #Set up a list to store the data
    if shapeType==0:
        timeSteps=230 #The length of time before mineralisation you care about
    else:
        timeSteps=1
    
    noVariables=18 #The number of variables you save
    
    input_rotation_filename = "Muller_gplates/Global_EarthByte_230-0Ma_GK07_AREPS.rot"
    # input_rotation_filename = "/Users/nbutter/Geodata/PlateModel/Shephard_etal_ESR2013_Global_EarthByte_2013.rot"
    #Create a rotation model to rotate the points back in time
    file_registry = pygplates.FeatureCollectionFileFormatRegistry()
    rotation_feature_collection = file_registry.read(input_rotation_filename)
    rotation_model = pygplates.RotationModel([rotation_feature_collection])
    
    #Initialise the array to store the data
    coregData=numpy.zeros((Nshp,timeSteps,noVariables))
    
    for i, nshpSub in enumerate(xrange(Nshp)):
            #!!! Warning: These are shapefile specific attributes, 
            #if you use a different shapefile, you must change these hard codes
            #to the corresponding age and plateID column in the shapefile
            if shapeType==0:
                shapeArray=numpy.array(shapes[nshpSub])
                age=0
                plateID=201
            else:
                shapeArray=numpy.array(shapes[nshpSub].points)
            #Here age is in the last column, and plateID in 7th.

            #1, 2 For porcu/main_edited.shp use 43 for the plateID
            #1, 2 For porcu/main_edited.shp use 9 for the age and 42 for the random age
            #3, 4 For XYBer14_t2_ANDES.shp use 7 for plateID
            #3, 4 For XYBer14_t2_ANDES.shp use 6 for age and -1 for Random
            if shapeType==1:
                age=round(recs[nshpSub][9]) 
                plateID=recs[nshpSub][43] 
            elif shapeType==2:
                age=round(recs[nshpSub][42]) 
                plateID=recs[nshpSub][43] 
            elif shapeType==3:
                age=int(round(recs[nshpSub][6]))
                plateID=recs[nshpSub][7] 
            elif shapeType==4:
                age=int(round(recs[nshpSub][-1]))
                plateID=recs[nshpSub][7]

            print(i, age, plateID)

            #Loop through each time step in the plate model
            for time in xrange(0,230,1):
            
                #If the point was formed at the current time (or 10Myr prior) then we 
                #want to find the surrounding plate properties.
                #if  (time > (age-10)) and (time < (age+10)) or (time > (age+20)) and (time < (age+30)):
                if  (time > (age-1)) and (time < (age+timeSteps)):
                    t=time-age
                    # print t, time
                    
                    #A raster file, to coregister with, these points have already been rotated
                    rasterfile="Muller_etal_2016_AREPS_Agegrids_v1.11/netCDF_0-230Ma/EarthByte_AREPS_v1.11_Muller_etal_2016_AgeGrid-"+str(time)+".nc"
                    # rasterfile="/Users/nbutter/Geodata/AgeGrid_20111110_Seton_etal_ESR/agegrid_final_mask_"+str(time)+".grd"
                    [x,y,z]=gridRead(rasterfile)
        
                    #A vector file to coregister with, these points have already been rotated
                    f = readCSV("Muller_convergence/subStats_"+str(time)+".csv")   
                    # f = readCSV("/Users/nbutter/Geodata/CONVERGENCE/Shep07/subStats_"+str(time)+".csv")   
                    lonlat=f[:,0:2]
                    
        
                    #-------------------#
                    #Reconstruct a point to its birth position
                    #-------------------#
                    if shapeType==0:
                        latlonPoint = pygplates.LatLonPoint(shapeArray[1],shapeArray[0])
                    else:
                        latlonPoint = pygplates.LatLonPoint(shapeArray[0][1],shapeArray[0][0])
                    point_to_rotate = pygplates.convert_lat_lon_point_to_point_on_sphere(latlonPoint)
            
                    finite_rotation = rotation_model.get_rotation(time, plateID)
                    birthPosition = finite_rotation * point_to_rotate
            
                    latlonBirth = pygplates.convert_point_on_sphere_to_lat_lon_point(birthPosition)
                    #allPoints = finite_rotation * pointSet_to_rotate
            
                    latBirth = latlonBirth.get_latitude()
                    lonBirth = latlonBirth.get_longitude()
        
                    #-------------------#
                    #Set the points for coregistering
                    region=5.0 #degrees
                    
                    #-------------------#
                    #Coregisterring raster 1
                    #-------------------#
                    #Find the region in index units
                    r=numpy.round(region/(x[1]-x[0]))
        
                    #Find the index unit of lat and lon
                    idxLon = (numpy.abs(x-lonBirth)).argmin()
                    idxLat = (numpy.abs(y-latBirth)).argmin()
                    
                    #Raster 1
                    c2=coregRaster([idxLon,idxLat],z,r)
                    #Hack to search further around the age grid if it can't find a match, note index units (not degrees)
                    # if numpy.isnan(c2):
                    #     print "Trying raster region: ", r+100.0
                    #     c2=coregRaster([idxLon,idxLat],z,r+100.0)
                    if numpy.isnan(c2):
                        c2=coregRaster([idxLon,idxLat],z,r+150.0)
                        print("Trying raster region: ", r+150.0)

                    #-------------------#
                    #Coregisterring vector 1
                    #-------------------#
                    index=coregPoint([lonBirth,latBirth],lonlat,region)
                    # if index=='inf':
                    #                         print "trying index region", region+10
                    #                         index=coregPoint([lonBirth,latBirth],lonlat,region+10.0)
                    if index=='inf':
                        print("trying index region", region+15)
                        index=coregPoint([lonBirth,latBirth],lonlat,region+15.0)
        
                    if numpy.isnan(c2) or index=='inf':
                        #if we have some null data, let's save it anyway, see what happens
                        #allData.append([shapeArray[0][0],shapeArray[0][1],0,0,\
                        #    age,age-time,0,0,0,0])
                        print("Skipping:", i, age, t, time, c2, index, lonBirth,latBirth)
        
                    else:
                        #Vector 1
                        segmentLength=f[index,3]
                        slabLength=f[index,9]
                        distSlabEdge=f[index,15]
                        
                        SPcoregNor=f[index,4]
                        SPcoregPar=f[index,12]
                        OPcoregNor=f[index,5]
                        OPcoregPar=f[index,13]
                        CONVcoregNor=f[index,10]
                        CONVcoregPar=f[index,11]
                        
                        subPolCoreg=f[index,8]
                        subOblCoreg=f[index,7]
                        
                        
#                         0 longitude, 
#                         1 latitude, 
#                         2 convRate(mm/yr), 
#                         3 segmentLength(km), 
#                         4 subPlateVel(mm/yr), 
#                         5 opVelocity(mm/yr), 
#                         6 trenchVel(mm/yr), 
#                         7 subObliquity(degrees), 
#                         8 subPolarity(degrees), 
#                         9 slabDistance(km), 
#                         10 OrthogonalConvergence(mm/yr),
#                         11 ParallelConvergence(mm/yr),
#                         12 ParallelSubPlateVel(mm/yr),
#                         13 ParallelOPvel(mm/yr),
#                         14 ParallelTrenchVel(mm/yr), 
#                         15 DistanceToSlabEdge(km),
#                         16 SPid, 
#                         17 TRENCHid, 
#                         18 OPid

                        if shapeType==0:
                            coregData[i,t,:]=[shapeArray[0],shapeArray[1],lonBirth,latBirth,\
                                        age,t,c2,\
                                        segmentLength,slabLength,distSlabEdge,\
                                        SPcoregNor,SPcoregPar,OPcoregNor,OPcoregPar,CONVcoregNor,CONVcoregPar,\
                                        subPolCoreg,subOblCoreg]

                        else:
                            coregData[i,t,:]=[shapeArray[0][0],shapeArray[0][1],lonBirth,latBirth,\
                                        age,t,c2,\
                                        segmentLength,slabLength,distSlabEdge,\
                                        SPcoregNor,SPcoregPar,OPcoregNor,OPcoregPar,CONVcoregNor,CONVcoregPar,\
                                        subPolCoreg,subOblCoreg]

    return(coregData)


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


def coregLoopTimeFirst(recs,shapes,fields,Nshp,shapeType):
    '''
    coregLoop

    This script reconstructs a shapefile of points to their birth time and 
    coregisters each point with another set of points, or a raster file.

    INPUTS:
    Can best be used with the output of readTopologyPlatepolygonFile() 
    which reads a shapefile and returns the 'records', 'shapes', 'fields',
    and 'number of shapes'.

    Hardcoded filenames/variables must be changed below:
        Rotation file - input_rotation_filename
        Time dependent raster files - rasterfile
        Time dependent vector csv outputs from 'convergence.py' - f

    OUTPUTS:
    Coregistered array: List of lat/lons with properties.

    METHOD:
    Takes a set of points (from a shapefile)
    Rotates the points back to their birth position
    Determines the point's birth position geophysical properties (coregisters)
    '''
    #Set up a list to store the data
    if shapeType==0:
        timeSteps=230 #The length of time before mineralisation you care about
    else:
        timeSteps=1

    noVariables=18 #The number of variables you save

    input_rotation_filename = "Muller_gplates/Global_EarthByte_230-0Ma_GK07_AREPS.rot"
    # input_rotation_filename = "/Users/nbutter/Geodata/PlateModel/Shephard_etal_ESR2013_Global_EarthByte_2013.rot"
    #Create a rotation model to rotate the points back in time
    file_registry = pygplates.FeatureCollectionFileFormatRegistry()
    rotation_feature_collection = file_registry.read(input_rotation_filename)
    rotation_model = pygplates.RotationModel([rotation_feature_collection])

    #Initialise the array to store the data
    coregData=numpy.zeros((Nshp,timeSteps,noVariables))

        #Loop through each time step in the plate model
    for time in xrange(0,230,1):

        #A raster file, to coregister with, these points have already been rotated
        rasterfile="Muller_etal_2016_AREPS_Agegrids_v1.11/netCDF_0-230Ma/EarthByte_AREPS_v1.11_Muller_etal_2016_AgeGrid-"+str(time)+".nc"
        # rasterfile="/Users/nbutter/Geodata/AgeGrid_20111110_Seton_etal_ESR/agegrid_final_mask_"+str(time)+".grd"
        [x,y,z]=gridRead(rasterfile)

        #A vector file to coregister with, these points have already been rotated
        f = readCSV("Muller_convergence/subStats_"+str(time)+".csv")   
        # f = readCSV("/Users/nbutter/Geodata/CONVERGENCE/Shep07/subStats_"+str(time)+".csv")

        lonlat=f[:,0:2]


        for i, nshpSub in enumerate(xrange(Nshp)):
                #!!! Warning: These are shapefile specific attributes, 
                #if you use a different shapefile, you must change these hard codes
                #to the corresponding age and plateID column in the shapefile
                if shapeType==0:
                    shapeArray=numpy.array(shapes[nshpSub])
                    age=0
                    plateID=201
                else:
                    shapeArray=numpy.array(shapes[nshpSub].points)
                #Here age is in the last column, and plateID in 7th.

                #1, 2 For porcu/main_edited.shp use 43 for the plateID
                #1, 2 For porcu/main_edited.shp use 9 for the age and 42 for the random age
                #3, 4 For XYBer14_t2_ANDES.shp use 7 for plateID
                #3, 4 For XYBer14_t2_ANDES.shp use 6 for age and -1 for Random
                if shapeType==1:
                    age=round(recs[nshpSub][9]) 
                    plateID=recs[nshpSub][43] 
                elif shapeType==2:
                    age=round(recs[nshpSub][42]) 
                    plateID=recs[nshpSub][43] 
                elif shapeType==3:
                    age=round(recs[nshpSub][6]) 
                    plateID=recs[nshpSub][7] 
                elif shapeType==4:
                    age=round(recs[nshpSub][-1]) 
                    plateID=recs[nshpSub][7]

                t=time-age
                print(i, age, time, t, plateID)

                #-------------------#
                #Reconstruct a point to its birth position
                #-------------------#
                if shapeType==0:
                    latlonPoint = pygplates.LatLonPoint(shapeArray[1],shapeArray[0])
                else:
                    latlonPoint = pygplates.LatLonPoint(shapeArray[0][1],shapeArray[0][0])
                
                point_to_rotate = pygplates.convert_lat_lon_point_to_point_on_sphere(latlonPoint)

                finite_rotation = rotation_model.get_rotation(time, plateID)
                birthPosition = finite_rotation * point_to_rotate

                latlonBirth = pygplates.convert_point_on_sphere_to_lat_lon_point(birthPosition)
                #allPoints = finite_rotation * pointSet_to_rotate

                latBirth = latlonBirth.get_latitude()
                lonBirth = latlonBirth.get_longitude()

                #-------------------#
                #Set the points for coregistering
                region=5.0 #degrees

                #-------------------#
                #Coregisterring raster 1
                #-------------------#
                #Find the region in index units
                r=numpy.round(region/(x[1]-x[0]))

                #Find the index unit of lat and lon
                idxLon = (numpy.abs(x-lonBirth)).argmin()
                idxLat = (numpy.abs(y-latBirth)).argmin()

                #Raster 1
                c2=coregRaster([idxLon,idxLat],z,r)
                #Hack to search further around the age grid if it can't find a match, note index units (not degrees)
                # if numpy.isnan(c2):
                #     print "Trying raster region: ", r+100.0
                #     c2=coregRaster([idxLon,idxLat],z,r+100.0)
                if numpy.isnan(c2):
                    c2=coregRaster([idxLon,idxLat],z,r+150.0)
                    print("Trying raster region: ", r+150.0)

                #-------------------#
                #Coregisterring vector 1
                #-------------------#
                index=coregPoint([lonBirth,latBirth],lonlat,region)
                # if index=='inf':
                #                         print "trying index region", region+10
                #                         index=coregPoint([lonBirth,latBirth],lonlat,region+10.0)
                if index=='inf':
                    print("trying index region", region+15)
                    index=coregPoint([lonBirth,latBirth],lonlat,region+15.0)

                if numpy.isnan(c2) or index=='inf':
                    #if we have some null data, let's save it anyway, see what happens
                    #allData.append([shapeArray[0][0],shapeArray[0][1],0,0,\
                    #    age,age-time,0,0,0,0])
                    print("Skipping:", i, age, t, time, c2, index, lonBirth,latBirth)

                else:
                    #Vector 1
                    segmentLength=f[index,3]
                    slabLength=f[index,9]
                    distSlabEdge=f[index,15]

                    SPcoregNor=f[index,4]
                    SPcoregPar=f[index,12]
                    OPcoregNor=f[index,5]
                    OPcoregPar=f[index,13]
                    CONVcoregNor=f[index,10]
                    CONVcoregPar=f[index,11]

                    subPolCoreg=f[index,8]
                    subOblCoreg=f[index,7]


    #                         0 longitude, 
    #                         1 latitude, 
    #                         2 convRate(mm/yr), 
    #                         3 segmentLength(km), 
    #                         4 subPlateVel(mm/yr), 
    #                         5 opVelocity(mm/yr), 
    #                         6 trenchVel(mm/yr), 
    #                         7 subObliquity(degrees), 
    #                         8 subPolarity(degrees), 
    #                         9 slabDistance(km), 
    #                         10 OrthogonalConvergence(mm/yr),
    #                         11 ParallelConvergence(mm/yr),
    #                         12 ParallelSubPlateVel(mm/yr),
    #                         13 ParallelOPvel(mm/yr),
    #                         14 ParallelTrenchVel(mm/yr), 
    #                         15 DistanceToSlabEdge(km),
    #                         16 SPid, 
    #                         17 TRENCHid, 
    #                         18 OPid

                    if shapeType==0:
                        coregData[i,t,:]=[shapeArray[0],shapeArray[1],lonBirth,latBirth,\
                                    age,t,c2,\
                                    segmentLength,slabLength,distSlabEdge,\
                                    SPcoregNor,SPcoregPar,OPcoregNor,OPcoregPar,CONVcoregNor,CONVcoregPar,\
                                    subPolCoreg,subOblCoreg]

                    else:
                        coregData[i,t,:]=[shapeArray[0][0],shapeArray[0][1],lonBirth,latBirth,\
                                    age,t,c2,\
                                    segmentLength,slabLength,distSlabEdge,\
                                  SPcoregNor,SPcoregPar,OPcoregNor,OPcoregPar,CONVcoregNor,CONVcoregPar,\
                                    subPolCoreg,subOblCoreg]

    return(coregData)