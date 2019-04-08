#!/usr/bin/env python
# coding: utf-8


#Import a few different tools and libraries. See Utils_coreg for more detail.
from Utils_coreg import gridRead, coregRaster, coregPoint

#Import data processing tools
import pygplates
import numpy

#andeanOutRand=open(r'Muller_Bertrand_coregistered_random.pkl','wb')
# pickle.dump(andesRand,andeanOutRand)
# andeanOutRand.close()

#This step require lots of time, so I load the modules then run these lines in a python terminal.
# f = readCSV("Muller_convergence/subStats_0.csv") #can use subStats_5, 10, etc for different transects
# lonlat=f[:,[0,1,17]] #Only lons, lats, and plate ID
# lonlats2 = lonlat[(lonlat[:,2])==201] #Pull out the trench ID with 201
# lonlats3=lonlats2[379:] #The shapefile has duplicate, overlapping lines, so just pull out the entire SAM trench
# Nshp = numpy.size(lonlats3,0)
# #Run the coregistration
# a=coregLoopTimeFirst(0,lonlats3,0,Nshp,0)
# andesPresent=numpy.array(a)
# andeanOutSamp=open(r'Muller_Bertrand_coregistered_sampleMuller0.pkl','wb')
# pickle.dump(andesPresent,andeanOutSamp)
# andeanOutSamp.close()



def coregLoop(sampleData,agegrid,kinfolder,input_rotation_filename):
    '''
    coregLoop

    This script reconstructs a shapefile of points to their birth time and 
    coregisters each point with another set of points, or a raster file.

    METHOD:
    Takes a set of points 
    Rotates the points back to their birth position
    Determines the point's birth position geophysical properties (coregisters)
    '''
    #Set up a list to store the data

    timeSteps=1
    noVariables=18 #The number of variables you save
    Nshp=len(sampleData)

    #Initialise the array to store the data
    coregData=numpy.zeros((Nshp,timeSteps,noVariables))

    #Create a rotation model to rotate the points back in time
    file_registry = pygplates.FeatureCollectionFileFormatRegistry()
    rotation_feature_collection = file_registry.read(input_rotation_filename)
    rotation_model = pygplates.RotationModel([rotation_feature_collection])


    #Loop over all the sample, coregistering each one.
    for i,nshpSub in enumerate(sampleData):

        lat=nshpSub[0]
        lon=nshpSub[1]
        age=int(nshpSub[2])
        plateID=int(nshpSub[3])
        
        print("Deposit:", i, "of", Nshp, "Lat:", lat, "Lon:", lon, "Age:", age,"PlateID:",plateID)

        #Loop through each time step in the plate model
        for time in xrange(0,230,1):
        
            #If the point was formed at the current time (or 10Myr prior) then we 
            #want to find the surrounding plate properties.
            #if  (time > (age-10)) and (time < (age+10)) or (time > (age+20)) and (time < (age+30)):
            if  (time > (age-1)) and (time < (age+timeSteps)):
                t=time-age
                # print t, time
                
                #A raster file, to coregister with, these points have already been rotated
                rasterfile=agegrid+str(time)+".nc"
                [x,y,z]=gridRead(rasterfile)
    
                #A vector file to coregister with, these points have already been rotated
                
                f=numpy.loadtxt(kinfolder+"subStats_"+str(time)+".csv", delimiter=',')  
                lonlat=f[:,0:2]
                
    
                #-------------------#
                #Reconstruct a point to its birth position
                #-------------------#
                latlonPoint = pygplates.LatLonPoint(lat,lon)
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
                if numpy.isnan(c2):
                    c2=coregRaster([idxLon,idxLat],z,r+150.0)
                    print("Trying raster region: ", r+150.0)

                #-------------------#
                #Coregisterring vector 1
                #-------------------#
                index=coregPoint([lonBirth,latBirth],lonlat,region)
                if index=='inf':
                    print("trying index region", region+15)
                    index=coregPoint([lonBirth,latBirth],lonlat,region+15.0)
    
                if numpy.isnan(c2) or index=='inf':
                    print("Skipping:", nshpSub, age, t, time, c2, index, lonBirth,latBirth)
    
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
                    
                    coregData[i,t,:]=[lon,lat,lonBirth,latBirth,                                        age,t,c2,                                        segmentLength,slabLength,distSlabEdge,                                        SPcoregNor,SPcoregPar,OPcoregNor,OPcoregPar,CONVcoregNor,CONVcoregPar,                                        subPolCoreg,subOblCoreg]
    
    #Return the filled coregistered array
    return(coregData)








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
        [x,y,z]=gridRead(rasterfile)

        #A vector file to coregister with, these points have already been rotated
        f=numpy.loadtxt(kinfolder+"subStats_"+str(time)+".csv", delimiter=',')   

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
                if numpy.isnan(c2):
                    c2=coregRaster([idxLon,idxLat],z,r+150.0)
                    print("Trying raster region: ", r+150.0)

                #-------------------#
                #Coregisterring vector 1
                #-------------------#
                index=coregPoint([lonBirth,latBirth],lonlat,region)
                if index=='inf':
                    print("trying index region", region+15)
                    index=coregPoint([lonBirth,latBirth],lonlat,region+15.0)

                if numpy.isnan(c2) or index=='inf':
                    #if we have some null data, let's save it anyway, see what happens
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
