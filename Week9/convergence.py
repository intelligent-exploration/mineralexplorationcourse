#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
201402018 Nathaniel Butterworth
convergence.py (Using pygplates r14878/GPlates 1.3)
Updated March 2017, working with pygplates r14878/GPlates 2.0
Updated March 2019, working with pygplates beta-revision-12

Determines convergence rates at subduction zones.
Uses output from GPlates plate reconstruction software.
Available from www.gplates.org

INPUT:
    
    Rotation file, Subduction zone shapefile Left, 
    Subduction zone shapefile Right, Plate polygon shapefile, 
    time to reconstruct from (probably should be same as age as polygon files),
    time to reconstruct to (probably should be timeFrom-1.)
    
    Can be run as a script, e.g:

for age in {0..10}; do echo ${age}; 
python convergence.py Shephard_etal_ESR2013_Global_EarthByte_2013.rot 
topology_subduction_boundaries_sL_${age}.00Ma.shp 
topology_subduction_boundaries_sR_${age}.00Ma.shp 
topology_platepolygons_${age}.00Ma.shp ${age}; done
    
    Or use run_convergence.py
    

Output
    longitude, latitude, convRate(mm/yr), 
    segmentLength(km), subPlateVel(mm/yr), opVelocity(mm/yr), 
    trenchVel(mm/yr), subObliquity(degrees), subPolarity(degrees), 
    slabDistance(km), OrthogonalConvergence(mm/yr), ParallelConvergence(mm/yr),
    ParallelSubPlateVel(mm/yr),ParallelOPvel(mm/yr),ParallelTrenchVel(mm/yr), 
    DistanceToSlabEdge(km),SPid, TRENCHid, 
    OPid

Method:
*Loads data. 
*For each point along a subduction zone determines the overriding (conjugate) 
plate and the subducting plate (the plate attached to the slab).
*Interpolates the data so it is evenly spaced.
*Then determines the relative and absolute motions of each plate.
*Convergence velocity is the relative rotation between the trench and SP.
*Absolut velocity is the motion of the SP.

To develop plate kinematic data from GPlates you may use 'convergence.py' 
script after following these directives.
#Load in plate polygons and rotation file to GPlates. Export these data:
#1. Export Resolved Topologies (CitcomS specific)
#2. Shapefiles
#3. ONLY select "Export all Plate Polygons to a single file"
#and "Export Plate Polygon segments to files based on feature type..."

#Then you can run this python scriptto produce a series of csv files 
with convergence data. 



Warnings:
*This script uses the 'new' (> Gplates 1.2) shapefile exports -> 
Resolved Topologies (CitcomS specific). "Export all Plate Polygons to a single 
file" and "Export Plate Polygon segments to files based on feature type".
Do not split polygons at the dateline (if you do the interpolation breaks.)
The 'old' shapefiles contain fewer shapefile attributes so the ordering is 
different. This is hardcoded in the script, but can be changed.

LIST OF FUNCTIONS:

readTopologyPlatepolygonFile
Reads shapefiles and returns the all the data fields

latlon2pygplates
Convert lat lon to the pygplates formats

subLoop
Determines the overriding and subducting plate ID along a trench.
Uses all the subduction zone data and plate polygon data.

kinloop
Determines the convergence and absolute velocities.

convloop
Determines the convergence rate. Needs two sequential points (to determine 
orthogonal direction) and an XYZ velocity.

convert180to360

convert360to180

calculateDistance

interpolatePoints

resampleData

reformatData

TODO:
*Not use shapefile output, instead use pygplates to get SZs (LH, RH, ID) and 
polygons directly from GPML feature collections.
*Deal with Triple Junctions.
*Plot convergence through time - done! 'coregTime.py'


   Copyright (C) 2014 The University of Sydney, Australia

   This file is part of GPlates.

   GPlates is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License, version 2, as published by
   the Free Software Foundation.

   GPlates is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

'''

import pygplates
print "Imported pygplates."
import shapefile
print "Imported shapefile."
import numpy
print "Imported numpy."
import sys
from scipy.interpolate import splprep, splev
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap, cm

### ---------------- INPUTS ---------------- ###
#Resolution of interpolation (in degrees)
global interpolation
interpolation=0.2

#input_rotation_filename = "/Users/nbutter/Geodata/PlateModel/Shephard_etal_ESR2013_Global_EarthByte_2013.rot"
#shapeLeft="/Users/nbutter/Geodata/PlateModel/Shephard/topology_subduction_boundaries_sL_75.00Ma.shp"
#shapeRight="/Users/nbutter/Geodata/PlateModel/Shephard/topology_subduction_boundaries_sR_75.00Ma.shp"
#shapeTopo="/Users/nbutter/Geodata/PlateModel/Shephard/topology_platepolygons_75.00Ma.shp"
#timeFrom=76
#timeTo=75
#outputFile='dummy.txt'



input_rotation_filename=sys.argv[1]
print "Rot:", input_rotation_filename

shapeLeft=sys.argv[2]
print "Left: ", shapeLeft
shapeRight=sys.argv[3]
shapeTopo=sys.argv[4]
timeFrom=int(sys.argv[5])
print(timeFrom)
# interpolation=float(sys.argv[6])
#
##shapeTopo1=sys.argv[7]
#
##If at present day, we rotate to present day (and not -1)
timeTo=timeFrom-1
if timeFrom==0:
   timeTo=0
#
outputFile='subStats_'+str(timeFrom)+'.csv'




### ---------------- Utility functions ---------------- ###
def readTopologyPlatepolygonFile(filename,typeSLT):
    '''
    Reads shapefiles and returns the all the data fields
    '''
    shapeRead = shapefile.Reader(filename)

    recs    = shapeRead.records()
    shapes  = shapeRead.shapes()
    fields  = shapeRead.fields
    Nshp    = len(shapes)
    
    return(recs,shapes,fields,Nshp)

def latlon2pygplates(lat,lon):
    '''
    Convert lat lon to the pygplates formats
    '''
    pointLatLon = pygplates.LatLonPoint(lat,lon)
    pointXYZ = pygplates.convert_lat_lon_point_to_point_on_sphere(pointLatLon)
    pointXYZcart = numpy.array([pointXYZ.get_x(), pointXYZ.get_y(), pointXYZ.get_z()]) 
    
    return(pointXYZ,pointXYZcart)
    

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

def calculateDistance(arrayoflatlons):
    '''
    Adds up the length of a set of geographical points that describe a line
    
    INPUT:
        [[x0,x1,x2...],[y0,y1,y2...]]
    OUTPUT:
        length (float)
    '''

    #Loop through all the points to make a line feature
    points=[]
    for i in arrayoflatlons:
        pointXYZ,pointXYZcart = latlon2pygplates(i[1],i[0])
        points.append(pointXYZ)

    #Make the line and find its distance    
    polyline = pygplates.PolylineOnSphere(points)
    length=polyline.get_arc_length()

    return length
    

##################################################
## SET of functions to interpolate subLoop data ##
##################################################
def unique(a):
    #Returns unique pairs in a 2d array
    #http://stackoverflow.com/a/8564438
     #These 2 lines preserve order
     #order = numpy.lexsort(a.T)
     #a = a[order]
     difference = numpy.diff(a, axis=0)
     ui = numpy.ones(len(a), 'bool')
     ui[1:] = (difference != 0).any(axis=1) 
     return a[ui]

       
def interpolatePoints(xs,ys,length):
    '''
    Interpolates points along an x-y (or lat-lon) line into a regulary spaced
    array. Resamples the number of points based on length of line. 
    INPUT:
        xs - [long0,long1,long2....]
        ys - [lat0,lat1,lat2...]
        length - float
    OUTPUT:
        Resampled xs - [newlon0,newlon1....]
        Resampled ys - [newlat0,newlat1....]
        
    TODO:
        
    Warnings:
        Uses a global variable called 'interpolation', which is the size of the
        interpolation in degrees
    '''
    from scipy.interpolate import splprep, splev    
    
    #Resolution of interpolation (in degrees)
    if 'interpolation' in globals():
        interDist=interpolation
    else:
        print 'No global interpolation, setting to 1 degree'
        interDist=1

    #The interpolation does not do well if the interpolation is over the +-180
    #So check if we need to convert it to 0-360
    #If the difference in lon values is very large, I assume interpolation fails
    lonDiff=numpy.diff(xs)
    #Check if each difference is bigger than 180 degrees (arbitrary)
    if any(abs(x) > 180.0 for x in lonDiff):
        tempxs=convert180to360(xs)
    else:
        tempxs=xs

    #Plot the original data set
    #print max(abs(lonDiff))
    #plt.scatter(tempxs,ys,c='b')

    #Interpolate the points
    tckp,u = splprep([tempxs,ys],k=1,s=0)
    u = numpy.arange(0,1,interDist/length)
    [xnew,ynew] = splev(u,tckp)

    
    #Now if it was converted, change it back to -180 to 180
    if any(abs(x) > 180.0 for x in lonDiff):
        tempxnew=convert360to180(xnew)
    else:
        tempxnew=xnew
    
    #Plot the interpolated one
    #plt.scatter(tempxnew,ynew,s=20,c='r',marker='x')
    
    #TODO: if you need to interpolate data also
    #f = interpolate.interp2d(xs, ys, zs, kind='linear')
    #znew = f(xnew, ynew)

    #plt.show()

    return [tempxnew,ynew]
    #return 

#Not used currently, but could be if we wanted to interpolated the data
#def interpolateZvalues(interpedX,interpedY,z):
#
#    f = interpolate.interp2d(interpedX, interpedY, z, kind='linear')
#    znew = f(xnew, ynew)
#
#    return znew

def resampleData(segments,SPid,TRENCHid,OPid,reformatArray):
    #For each segment, interpolate the data
    segment=numpy.array(segments)

    if len(segment) > 1:
        #check the numbers in the array are unique
    	f=numpy.dstack((segment[:,0],segment[:,1]))
        g=numpy.reshape(f,(f.shape[1],2))
        h=unique(g)

        if len(h[:,0]) > 1:
            length=calculateDistance(h)
            reformatArray.append(['> Length:' , length*6371.0])
            [interpLon,interpLat]=interpolatePoints(h[:,0],h[:,1],numpy.rad2deg(length))
            oneSP=numpy.ones(len(interpLon))*SPid
            oneTRENCH=numpy.ones(len(interpLon))*TRENCHid
            oneOP=numpy.ones(len(interpLon))*OPid

            #Add the new resampled points to our resampled Array
            for i,k in enumerate(interpLon):
                reformatArray.append([interpLon[i],interpLat[i],oneSP[i],oneTRENCH[i],oneOP[i]])
        
    #Reset the segment count and go to the next set
    segments=[]
    return reformatArray


def reformatData(subArray):
    '''
    Calculates the length of subduction segment and resamples the data
    INPUT: Array from subloop of the form:
        [Lon, Lat, SPid, TRENCHid, OPid], with headers indicated by '>'
    OUTPUT:
        Resampled data of the same form as input:
        [Lon, Lat, SPid, TRENCHid, OPid], with headers indicated by '>'
    '''
    #print subArray
    #print 'The initial array is: ', len(subArray)
    reformatArray=[]
    start=1
    secondPoint=1 
    end=0   
    #Go through all the data
    segments=[]
    for index, point in enumerate(subArray[:]): 
        #Check if the data is header (at the start/end) of a segment
        if isinstance(point[0],str):
            secondPoint=1
            #Save the data as we go if it is a header
            if start:
                reformatArray.append(point)

            #now check if we need to interpolate the data
            #This should only be reached at the end-of-segment (not start)
            if (point[0].find('Name') != -1 and not start):
                reformatArray=resampleData(segments,SPid,TRENCHid,OPid,reformatArray)
                segments=[]
                
            #Append the 'header' information at the correct time.    
            if not start:
                reformatArray.append(point)

        #If we have data append it to the array...
        if isinstance(point[0],float):
            #Store the plate ids as we go, because we need them when we go back 
            #to the end-of-segment conditional if, above.
            SPid=subArray[index][2]
            TRENCHid=subArray[index][3] 
            OPid=subArray[index][4]
            
            #Check if the plate ids change, so we can resample that segment
            if secondPoint==0:
                SPid2=subArray[index-1][2]
                TRENCHid2=subArray[index-1][3] 
                OPid2=subArray[index-1][4]
                secondPoint=1

                if (SPid != SPid2 or TRENCHid != TRENCHid2 or OPid != OPid2):
                    reformatArray=resampleData(segments,SPid2,TRENCHid2,OPid2,reformatArray)
                    segments=[]
            
            #Add the point to the segment array, that we will use when we are 
            #back in the end-of-segment conditional if, above.
            segments.append(point)
        
            #If you are at the very last point in the set of points
            if index==len(subArray):
                reformatArray=resampleData(segments,SPid,TRENCHid,OPid,reformatArray)
            
            #Once you have a real point, you are waiting for the next header
            start=0
            secondPoint=0
                        
    return(reformatArray)
                


#########################################################
### ---------------- Major Functions ---------------- ###
#########################################################   
def subLoop(recsSub,shapesSub,fieldsSub,NshpSub,recsTopo,shapesTopo,fieldsTopo,NshpTopo,isright):
    '''
    Determines the overriding and subducting plate ID along a trench.
    Takes the subduction zone and plate polygon shapefile data as input.
    Returns a list of points of the form [Lon, Lat, SPid, TRENCHid, OPid]
    And puts some header information indicated by '>' for each subduction segment.
    '''
    #Initiialise the output array: [lon, lat, SPid, TRENCHid, OPid] 
    testingArray=[]
    FinalArray=[]
    
    #Go through subduction zones 1 by 1
    #print "Subdction zone number in file, and PlateID:"
    #NshpSub=[8,9]
    for nshpSub in xrange(NshpSub):
        
        #Get the Lat Lon points along the trench
        shapeArray=shapesSub[nshpSub].points
        
        #Reverse the order of LH arrays so convergence is +ve        
        reversedSub=0
        if not isright:
            #print "Reversing, LH array"
            shapeArray=shapeArray[::-1]
            reversedSub=1
        
        #Get the Trench ID and Time
        # Note this section looks both for a PLATE_ID and a PLATEID1 field for the trench
        # index, due to changes between different GPlates versions. If both are present for
        # some reason, the TRENCHid would be taken from the second one to be encountered
        for i in range(len(fieldsSub)):
            if 'PLATE_ID' in fieldsSub[i]:
                TRENCHidIndex = i-1
            elif 'PLATEID1' in fieldsSub[i]:
                TRENCHidIndex = i-1
            elif 'TIME' in fieldsSub[i]:
                TimeIndex = i-1
            elif 'FEATURE_ID' in fieldsSub[i]:
                GplatesTagIndex = i-1
            elif 'NAME' in fieldsSub[i]:
                NameIndex = i-1

        #Get the Trench ID and Time
        TRENCHid=recsSub[nshpSub][TRENCHidIndex] #3 new, 0 old (GPlates 1.2)        
        Time=recsSub[nshpSub][TimeIndex] #2 old
        GplatesTag=recsSub[nshpSub][GplatesTagIndex]
        name=recsSub[nshpSub][NameIndex]
        #print nshpSub,TRENCHid,GplatesTag,name


        #Store the Lat and Lon in seperate arrays, remove end points because 
        #they are prone to errors and confusing triple junctions!
        shapeArrayLon = [round(i[0],4) for i in shapeArray[:]]
        shapeArrayLat = [round(i[1],4) for i in shapeArray[:]]
        
        shapeArray=zip(shapeArrayLon,shapeArrayLat)

        #Append some header information
        #NPB!!! This is currently not used for much
        FinalArray.append(['> Name:' , name])
        FinalArray.append(['> GplatesTag:' , GplatesTag])
        FinalArray.append(['> PlateId:' , TRENCHid])
        FinalArray.append(['> Age:' , Time])
        FinalArray.append(['> Begin:' , 0])
        FinalArray.append(['> End:' , 0])
        FinalArray.append(['> Type:' , 0])
        FinalArray.append(['> Polygon:' , 0])
        FinalArray.append(['> Reversed:' , reversedSub])
        FinalArray.append(['> BoundaryPoints:' , len(shapeArray)])
        
        #Now we must loop through every point along the subduction zone and
        #check the polygons that it is attached to.
        for index, point in enumerate(shapeArray):
            OPcount=0
            SPcount=0
            for nshpTopo in xrange(NshpTopo):   
                #Get the Lat Lon points of the plate polygon
                topoArray=shapesTopo[nshpTopo].points

                #Store the Lat and Lon in seperate arrays
                topoArrayLon = [round(i[0],4) for i in topoArray]
                topoArrayLat = [round(i[1],4) for i in topoArray]

                topoArray=zip(topoArrayLon,topoArrayLat)
                
                #topoArray=numpy.around(shapesTopo[nshpTopo].points,decimals=5)
                
                #Test whether the subduction point is in the current
                #polygon. If it is, we can find out more information...
                try:
                    i = topoArray.index(point)
                    pointFlag = 1
                except:
                    pointFlag = 0
                
                #########################################
                #Now if the point exists on the polygon #
                #Let's find out whether it is OP or SP  #
                #########################################
                if pointFlag == 1:
                    #Get the PLATEid of the polygon (same comment as above - looks for field called either PLATE_ID or PLATEID1)
                    for i in range(len(fieldsTopo)):
                        if 'PLATE_ID' in fieldsTopo[i]:
                            PlateIndex = i-1
                        elif 'PLATEID1' in fieldsTopo[i]:
                            PlateIndex = i-1
                    PLATEid=recsTopo[nshpTopo][PlateIndex] #3 new, 0 old (GPlates 1.2) 

                    #Now get the entire polygon as a pygplates entitiy
                    latlonPoints=[]
                    for latlon in topoArray:
                        pointLatLon = pygplates.LatLonPoint(latlon[1],latlon[0])
                        latlonPoints.append(pygplates.convert_lat_lon_point_to_point_on_sphere(pointLatLon))

                    polygon = pygplates.PolygonOnSphere(latlonPoints)
                                        
                    #Convert our trench segment to pygplates format
                    if index < len(shapeArray)-1:
                        latPoint1 = shapeArray[index][1]
                        lonPoint1 = shapeArray[index][0]
                        pointXYZ, pointXYZcart = latlon2pygplates(latPoint1,lonPoint1)
                    
                        latPoint2 = shapeArray[index+1][1]
                        lonPoint2 = shapeArray[index+1][0]
                        pointXYZ2, pointXYZcart2 = latlon2pygplates(latPoint2,lonPoint2)
                    #If we are at the last trench segment, just jump back a point
                    else:
                        latPoint1 = shapeArray[index-1][1]
                        lonPoint1 = shapeArray[index-1][0]
                        pointXYZ, pointXYZcart = latlon2pygplates(latPoint1,lonPoint1)

                        latPoint2 = shapeArray[index][1]
                        lonPoint2 = shapeArray[index][0]
                        pointXYZ2, pointXYZcart2 = latlon2pygplates(latPoint2,lonPoint2)

                        
                    #Find the vector between 2 points on the trench
                    line = pygplates.PolylineOnSphere([pointXYZ,pointXYZ2])

                    #Find the midpoint between those two points
                    midPoint = line.get_centroid() 
                    midLatLon = pygplates.convert_point_on_sphere_to_lat_lon_point(midPoint)
                    midPointLatLon =  numpy.array([midLatLon.get_longitude(), midLatLon.get_latitude()])
                    
                    #Find the bearing of point 1 to point 2 (ALL IN RADIANS)
                    #eqns from http://www.movable-type.co.uk/scripts/latlong.html
                    radLat1 = numpy.radians(latPoint1)
                    radLat2 = numpy.radians(latPoint2)
                    radLon1 = numpy.radians(lonPoint1)
                    radLon2 = numpy.radians(lonPoint2)
                            
                    #θ = atan2( sin(Δλ).cos(φ2), cos(φ1).sin(φ2) − sin(φ1).cos(φ2).cos(Δλ) )
                    bearingRad = numpy.arctan2( numpy.sin(radLon2-radLon1) * numpy.cos(radLat2),\
                        numpy.cos(radLat1)*numpy.sin(radLat2) - numpy.sin(radLat1)*numpy.cos(radLat2)*numpy.cos(radLon2-radLon1))
                    #Now we march inside the polygon, to check if it is OP or SP
                    #Turn 90 degrees to the left or right
                    #if isright:
                    bearingPerp = bearingRad + 1.57
                    #else:
                    #    bearingPerp = bearingRad - 1.57
                    
                    #March into the plate a few km along the bearing.
                    #new Lat Lon are given below
                    d = 20.0 #km #Arbitrary small distance...
                    R=6371.0 #km Radius of Earth km
                    
                    #φ2 = asin( sin(φ1)*cos(d/R) + cos(φ1)*sin(d/R)*cos(θ) )
                    lat2 = numpy.arcsin(numpy.sin(radLat1)*numpy.cos(d/R) + numpy.cos(radLat1)*numpy.sin(d/R)*numpy.cos(bearingPerp))
                    #λ2 = λ1 + atan2( sin(θ)*sin(d/R)*cos(φ1), cos(d/R)−sin(φ1)*sin(φ2) )
                    lon2 = radLon1 + numpy.arctan2(numpy.sin(bearingPerp)*numpy.sin(d/R)*numpy.cos(radLat1),numpy.cos(d/R)-numpy.sin(radLat1)*numpy.sin(radLat2))        
                    pointSurface = pygplates.LatLonPoint(numpy.degrees(lat2),numpy.degrees(lon2))
                    pointSurfaceXYZ = pygplates.convert_lat_lon_point_to_point_on_sphere(pointSurface)
                    
                    #Test whether it is in the polygon or not
                    cwcheck = polygon.is_point_in_polygon(pointSurfaceXYZ)

                    #If it is then we have marched into the cojugate plate (probably)!
                    if cwcheck:
                        OPid=PLATEid
                        OPcount+=1
                    else:
                        SPid=PLATEid
                        SPcount+=1
                        
                #If the point does not exist on the polygon, pass                    
                else:
                    pass
                
            #After looping through the polygons, we should have one OPid, 
            #and one SPid, but these may be wrong for triple junctions.            
            if (OPcount == 1) and (SPcount == 1):
                #Now exit the Polygon loop and save the information about the point and move on 
                #print point[0],point[1], SPid, TRENCHid, OPid
                FinalArray.append([point[0],point[1], SPid, TRENCHid, OPid])
            else:
                pass
                #print OPcount, SPcount
                #print "Warning, point: ", point, TRENCHid, " found no match in polygons."
                #print "'BoundaryPoints' will not match actual points!"
    #print FinalArray            
    return(FinalArray)


def kinloop(timeFrom,timeTo,data,input_rotation_filename):
    '''
    Determines the absolute and relative velocities of the plates and trench.

    INPUT (with examples):
        
    timeFrom=2 #integer
    timeTo=1 #integer
    data = [[14.7703, 43.5369, 701, 301, 301], [18.8934, 41.5683, 701, 301, 301], 
            [18.6062, 39.7998, 701, 301, 301], [19.4349, 39.3242, 701, 301, 301]]
            #List of points [[Lon, Lat, SPid, TRENCHid, OPid]]
    input_rotation_filename='Shephard_etal_ESR2013_Global_EarthByte_2013.rot' #String
    
    OUTPUT: (array equivalent to 'data' plus Absoulte XYZ and Convergence XYZ 
        velocity components.)

    [14.7703, 43.5369, 701, 301, 301, -13.02077,10.4195,10.45488,-4.86172,
        -4.62026,6.18716]
        
    WARNING: If this function is used stand-alone, without an idea of left
        or right, then it may be backward. Points should be entered CW.

    '''
    #Set the times of interest
    #timeFrom=3 #from xx Ma
    #timeTo=2 # to xx Ma
    
    #Read the SubdctionZoneAttributes data
    kinematicsArray=data

    #Read the rotation file and build a rotation model    
    #input_rotation_filename = "/Users/nbutter/DataMining/CODE/SubductionZoneAttributes/oldTestdata/Global_EarthByte_TPW_CK95G94_PP_20110412.rot"
    file_registry = pygplates.FeatureCollectionFileFormatRegistry()
    rotation_feature_collection = file_registry.read(input_rotation_filename)
    rotation_model = pygplates.RotationModel([rotation_feature_collection])

    #Build the rotations for the times we want    
    fromTimeRotation=rotation_model.get_reconstruction_tree(timeFrom)
    toTimeRotation=rotation_model.get_reconstruction_tree(timeTo)
    
    #Make an array to store our data in
    speedArray=[]
    reversedSub=0 #Initilise this here in case the data does not have headers
    
    #Loop through all the SubdctionZoneAttributes data and find velocities
    for index, point in enumerate(kinematicsArray[:]): 
        #print point[0]
        #See if the point is not a header, or if the next point is not a header
        if isinstance(point[0],str):
            #Also, if it is a header, get the order of the points
            if point[0].find('Reversed') != -1:
                reversedSub=int(point[-1])
                #print reverse
        if len(point) == 5:
            lon = float(point[0])
            lat = float(point[1])
            OPid = int(point[4])
            TRENCHid = int(point[3])
            SPid = int(point[2])

            #Convert point LatLon to XYZ
            pointXYZ, pointXYZcart = latlon2pygplates(lat,lon)
            
            #Get the Euler poles for each of the relevant rotations
            stageTrench = pygplates.ReconstructionTree.get_equivalent_stage_rotation(fromTimeRotation,toTimeRotation,TRENCHid)
            poletTrench,angleTrench = stageTrench.get_euler_pole_and_angle()
            poleTrench = [poletTrench.get_x(), poletTrench.get_y(), poletTrench.get_z()]

            #Absolut velocity is the motion of the SP
            stageSP = pygplates.ReconstructionTree.get_equivalent_stage_rotation(fromTimeRotation,toTimeRotation,SPid)
            poletSP,angleSP = stageSP.get_euler_pole_and_angle()
            #If you want the latitude and longitude, rather than the cartesian:
            #poleSPLatLon = pygplates.convert_point_on_sphere_to_lat_lon_point(poleSP)
            #SPLatLon = [poleSPLatLon.get_latitude(), poleSPLatLon.get_longitude()] 
            poleSP = [poletSP.get_x(), poletSP.get_y(), poletSP.get_z()]

            #Convergence is the relative rotation between the trench and SP
            #if reversedSub==1:
            #convRotation = pygplates.get_relative_stage_rotation(fromTimeRotation,toTimeRotation,TRENCHid,SPid)
            #else:
            convRotation = pygplates.ReconstructionTree.get_relative_stage_rotation(fromTimeRotation,toTimeRotation,SPid,TRENCHid)
            poletCONV,angleCONV = convRotation.get_euler_pole_and_angle()
            poleCONV = [poletCONV.get_x(), poletCONV.get_y(), poletCONV.get_z()]

            #The OP motion
            stageOP = pygplates.ReconstructionTree.get_equivalent_stage_rotation(fromTimeRotation,toTimeRotation,OPid)
            poletOP,angleOP = stageOP.get_euler_pole_and_angle()
            poleOP = [poletOP.get_x(), poletOP.get_y(), poletOP.get_z()]
            
            #Angular velocity (Degrees / Myr) #but need it in radians.. not used
            omegaCONV = numpy.degrees(angleCONV)/(timeFrom-timeTo)
            omegaSP = numpy.degrees(angleSP)/(timeFrom-timeTo) 

            #Velocity vectors at each point in (Earth Radius Units) / Myr
            absoluteV = numpy.cross(poleSP,pointXYZcart) * angleSP/(timeFrom-timeTo) 
            convergenceV = numpy.cross(poleCONV,pointXYZcart) * angleCONV/(timeFrom-timeTo) 
            opV = numpy.cross(poleOP,pointXYZcart) * angleOP/(timeFrom-timeTo)
            trenchV = numpy.cross(poleTrench,pointXYZcart) * angleTrench/(timeFrom-timeTo)            

            #Velocivites vectors in mm/yr (km/Myr)
            Rearth=6371.0 # m 
            vABS = absoluteV * Rearth
            vCONV = convergenceV * Rearth
            vOP = opV * Rearth
            vTrench = trenchV * Rearth

            speedArray.append([lon,lat,SPid,TRENCHid,OPid,\
                round(vABS[0],5),round(vABS[1],5),round(vABS[2],5),\
                round(vCONV[0],5),round(vCONV[1],5),round(vCONV[2],5),\
                round(vOP[0],5),round(vOP[1],5),round(vOP[2],5),\
                round(vTrench[0],5),round(vTrench[1],5),round(vTrench[2],5),\
                int(SPid), int(TRENCHid), int(OPid)])
            #print [lon,lat,SPid,TRENCHid,OPid,round(vABS[0],5),round(vABS[1],5),round(vABS[2],5),round(vCONV[0],5),round(vCONV[1],5),round(vCONV[2],5)]
        #If it is a header, just print it back out.
        else:
            #pass
            speedArray.append(point)

    return(speedArray)



def kinstats(kinematicsArray):
    '''
    Determines the convergence and other details of a subducting point
    Requires at least 2 points per trench segment do calculate the orthogonal
    direction.
    
    INPUT (output from kinloop) or:
    [Lon, Lat, SPid, TRENCHid, OPid, AbsX, AbsY, AbsZ, ConvX, ConvY, ConvZ]

    OUTPUT:
    
    longitude, latitude, convRate(mm/yr), segmentLength(km), subPlateVel(mm/yr), 
    opVelocity(mm/yr), trenchVel(mm/yr), subObliquity(degrees), 
    subPolarity(degrees), CumSlabLength(km), OrthogonalConvergence(mm/yr),
    ParallelConvergence(mm/yr),ParallelSubPlateVel(mm/yr),
    ParallelOPvel(mm/yr),ParallelTrenchVel(mm/yr), DistanceToSlabEdge(km)
    '''
    #Read the SubdctionZoneAttributes data
    #subData="/Users/nbutter/DataMining/CODE/SubductionZoneAttributes/oldTestdata/2.vel"
    #print "Reading ", subData  
    #kinematicsArray=readSZAoutput(subData)

    #Initiliase an output array
    convArray=[]
    
    #Set these variables here in case the data is not in the headers
    distEdge = 0 #Set the distance of the point to the end of the segment to 0
    slabLength = 0
    distEdgeTotal = 0
    for index, point in enumerate(kinematicsArray[:-1]): 
        #See if the point is not a header, or if the next point is not a header
        if isinstance(point[0],str):
            if point[0].find('Length') != -1:
                slabLength=float(point[-1])
                distEdgeTotal=0
                #print index

        if (len(point) > 2) & (len(kinematicsArray[index+1]) > 2):
            #Save the values we will need to use
            lon = float(point[0])
            lat = float(point[1])
            conv = numpy.array([float(point[8]),float(point[9]),float(point[10])])
          
            lon2 = float(kinematicsArray[index+1][0])
            lat2 = float(kinematicsArray[index+1][1])
            conv2 = numpy.array([float(kinematicsArray[index+1][8]),\
                float(kinematicsArray[index+1][9]),float(kinematicsArray[index+1][10])])

            vabs = numpy.array([float(point[5]),float(point[6]),float(point[7])])
            vabs2 = numpy.array([float(kinematicsArray[index+1][5]),\
                float(kinematicsArray[index+1][6]),float(kinematicsArray[index+1][7])])
            
            vop = numpy.array([float(point[11]),float(point[12]),float(point[13])])
            vop2 = numpy.array([float(kinematicsArray[index+1][11]),\
                float(kinematicsArray[index+1][12]),float(kinematicsArray[index+1][13])])

            vtrench = numpy.array([float(point[14]),float(point[15]),float(point[16])])
            vtrench2 = numpy.array([float(kinematicsArray[index+1][14]),\
                float(kinematicsArray[index+1][15]),float(kinematicsArray[index+1][16])])

            SPid = kinematicsArray[index+1][17]
            TRENCHid = kinematicsArray[index+1][18] 
            OPid = kinematicsArray[index+1][19]

            #Convert them to the pygplates formats we will need
            pointXYZ, pointXYZcart = latlon2pygplates(lat,lon)
            pointXYZ2, pointXYZcart2 = latlon2pygplates(lat2,lon2)

            #Find the difference length of the two points
            diffVec = pointXYZcart2 - pointXYZcart

            #Find the normal direction between the two points
            normal = numpy.cross(diffVec,pointXYZcart)
            normalise = numpy.linalg.norm(normal)
            normal2 = normal/normalise
            
            diffnorm = numpy.linalg.norm(diffVec)
            diffnorm2 = diffVec / diffnorm

            #Get the convergence rate of the segment
            avgConv = (conv + conv2) / 2.0
            #Conv rate was orthognal to the trench, and we took the abs value to 
            #account for places where the conv was negative, but now we use this
            #column to be the magnitude of conv rate.
            convRate = numpy.sqrt(numpy.dot(avgConv,avgConv)) #Orthogonal to the trench
            convPerp = numpy.dot(normal2,avgConv)
            convPar = numpy.dot(diffnorm2,avgConv) #Parallel to the trench

            avgAbs = (vabs + vabs2) / 2.0
            orthAbs = numpy.dot(normal2,avgAbs)
            parAbs = numpy.dot(diffnorm2,avgAbs)
            
            avgOP = (vop + vop2) / 2.0
            orthOP = numpy.dot(normal2,avgOP)
            parOP = numpy.dot(diffnorm2,avgOP)
            
            avgTrench = (vtrench + vtrench2) / 2.0
            orthTrench = numpy.dot(normal2,avgTrench)
            parTrench = numpy.dot(diffnorm2,avgTrench)
        
            
            #a = numpy.arccos(numpy.dot(a1,a2) / numpy.dot(numpy.abs(a1),numpy.abs(a2)))
            #Old method
            #subObliquity = numpy.arccos(numpy.dot(diffVec,avgConv) / numpy.dot(numpy.abs(diffVec),numpy.abs(avgConv)))
            #subObliquity = numpy.degrees(subObliquity)
            #New, better definition
            subObliquity = numpy.arccos(numpy.dot(diffVec,avgConv) / (numpy.linalg.norm(diffVec) * numpy.linalg.norm(avgConv)))
            #NPB: I don't think this should be absolute value ,as you lose information!
            #buuut, you can get it back by doing arctan(Parallel/Perpendicular)
            # Or you can out put it here with: subObliquitySigned = numpy.degrees(subObliquity) - 90)
            subObliquity = numpy.abs(numpy.degrees(subObliquity) - 90)
            

            #Find the bearing of point 1 to point 2 (ALL IN RADIANS)
            #eqns from http://www.movable-type.co.uk/scripts/latlong.html
            radLat1 = numpy.radians(lat)
            radLat2 = numpy.radians(lat2)
            radLon1 = numpy.radians(lon)
            radLon2 = numpy.radians(lon2)
                            
            #θ = atan2( sin(Δλ).cos(φ2), cos(φ1).sin(φ2) − sin(φ1).cos(φ2).cos(Δλ) )
            bearingRad = numpy.arctan2( numpy.sin(radLon2-radLon1) * numpy.cos(radLat2),\
                numpy.cos(radLat1)*numpy.sin(radLat2) - numpy.sin(radLat1)*numpy.cos(radLat2)*numpy.cos(radLon2-radLon1))

            #subObliquity = numpy.degrees(bearingRad)
            #Angle of the trench relative to North
            #NS and subducting toward the east is 0
            #NS and subducting west is 180/-180
            #i.e subduction toward the N=-90, S90, W-180/+180, E0
            subPolarity = numpy.degrees(bearingRad)

            #if subPolarity < 0:
            #    subPolarity = subPolarity+360.0
            
            #Get the length of segment between the two points
            Rearth=6371.0 # Radius of Earth in km
            line = pygplates.PolylineOnSphere([pointXYZ,pointXYZ2])
            distance = line.get_arc_length() * Rearth
            
            #If we are still in the same segment we can update the distance to the edge.            
            distEdge+=distance
            if distEdge < slabLength/2.0:
                distEdgeTotal+=distance
            else:
                distEdgeTotal-=distance
            #print slabLength/2.0, distance, distEdge, distEdgeTotal

            #Find the midpoint of line and convert to lat lon
            midPoint = line.get_centroid() 
            LatLonMidPoint = pygplates.convert_point_on_sphere_to_lat_lon_point(midPoint)
            lonMid, latMid =  LatLonMidPoint.get_longitude(), LatLonMidPoint.get_latitude()

            #Return the array, Longitude, Latitude, ConvergenceRate(mm/yr),distance(km)
            convArray.append([lonMid, latMid, convRate, distance,\
            orthAbs, orthOP, orthTrench,subObliquity,subPolarity,distEdge,\
            convPerp,convPar,parAbs,parOP,parTrench,distEdgeTotal,\
            int(SPid), int(TRENCHid), int(OPid)])

        #If it is a header, just print it back out.
        else:
            distEdge = 0 #Reset the segment distance to 0
            #pass #convArray.append(point)
            
    return(convArray)    

### ---------------- Main function ---------------- ###
def main(shapeTopo,shapeLeft,shapeRight,input_rotation_filename,timeFrom,timeTo,outputFile,\
        ):

    #Read in the shapefiles
    [recsSubL,shapesSubL,fieldsSubL,NshpSubL]=readTopologyPlatepolygonFile(shapeLeft,'L')
    [recsSubR,shapesSubR,fieldsSubR,NshpSubR]=readTopologyPlatepolygonFile(shapeRight,'R')
    [recsTopo,shapesTopo,fieldsTopo,NshpTopo]=readTopologyPlatepolygonFile(shapeTopo,'T')

    #Deforming plate polygons need extra polygon file to be read
    #try:
    #    [recsTopo1,shapesTopo1,fieldsTopo1,NshpTopo1]=readTopologyPlatepolygonFile(shapeTopo1,'T')  
    #    recsTopo=numpy.r_[recsTopo,recsTopo1] 
    #    shapesTopo=numpy.r_[shapesTopo,shapesTopo1]
    #    NshpTopo = NshpTopo+NshpTopo1
    #except:
    #    pass

    print "Reading Polygon files..."
    print shapeTopo
    
    print "Looping over LH subduction zones..."
    print shapeLeft
    subArrayL=subLoop(recsSubL,shapesSubL,fieldsSubL,NshpSubL,recsTopo,shapesTopo,fieldsTopo,NshpTopo,False)
    #print subArrayL
    #Do the same for the right hand
    print "Looping over RH subduction zones..."
    print shapeRight
    subArrayR=subLoop(recsSubR,shapesSubR,fieldsSubR,NshpSubR,recsTopo,shapesTopo,fieldsTopo,NshpTopo,True)
    
    #Append the two arrays
    subArray=subArrayL+subArrayR

    #Resample the data so the points are evenly interpolated
    if interpolation != 0:
        print 'Interpolating data to ',interpolation,' degrees...'
        reformArray = reformatData(subArray)
    
    #print reformArray
    #Still works if you do not resample the data!
    else:
        print 'Skipping interpolation.'
        reformArray = subArray
    
    #print reformArray
    #Determine the kinemastics (absoulte and relative velocities)
    print "Finding velocities of all the points..."
    kinArray = kinloop(timeFrom,timeTo,reformArray,input_rotation_filename)
    
    #Determine the convergence rates
    print "Finding trench orthogonal convergence..."
    convArray = kinstats(kinArray)
    
    #Return the convergence rate data
    return(convArray)


if __name__ == "__main__":
    
    print 'Running main...'
    FinalArray=main(shapeTopo,shapeLeft,shapeRight,input_rotation_filename,timeFrom,timeTo,outputFile)
    
    f=numpy.array(FinalArray)

    
    #Save the data out
    numpy.savetxt(outputFile, f, delimiter=',',fmt='%1.5f')
    print "Done ", outputFile
    
    #you can plot the data at this point like this
#    import matplotlib.pyplot as plt
#    from mpl_toolkits.basemap import Basemap, cm
#
#    pmap = Basemap(projection='hammer', lat_0=0, lon_0=180,
#               resolution='l')
#    pmap.drawmapboundary(fill_color='white')
#    pmap.fillcontinents(color='grey', lake_color='white', zorder=0)
#    pmap.drawmeridians(numpy.arange(0, 360, 30))
#    pmap.drawparallels(numpy.arange(-90, 90, 30))
#
#    xh, yh = pmap(f[:,0], f[:,1])
#    l1 = pmap.scatter(xh, yh, 15, marker='o', c=f[:,2],
#                  edgecolor='none', vmin=-100, vmax=100)
#    pmap.colorbar(l1,location='bottom',pad="5%")
#    plt.title('Convergence Rate (km/Myr)')
#    plt.show()
#
#    pmap = Basemap(projection='hammer', lat_0=0, lon_0=180,
#               resolution='l')
#    pmap.drawmapboundary(fill_color='white')
#    pmap.fillcontinents(color='grey', lake_color='white', zorder=0)
#    pmap.drawmeridians(numpy.arange(0, 360, 30))
#    pmap.drawparallels(numpy.arange(-90, 90, 30))
#
#    xh, yh = pmap(f[:,0], f[:,1])
#    l1 = pmap.scatter(xh, yh, 15, marker='o', c=f[:,10],
#                  edgecolor='none', vmin=-100, vmax=100)
#    pmap.colorbar(l1,location='bottom',pad="5%")
#    plt.title('Convergence Velocity (km/Myr)')
#    plt.show()
#    
#    
#    pmap = Basemap(projection='hammer', lat_0=0, lon_0=180,
#               resolution='l')
#    pmap.drawmapboundary(fill_color='white')
#    pmap.fillcontinents(color='grey', lake_color='white', zorder=0)
#    pmap.drawmeridians(numpy.arange(0, 360, 30))
#    pmap.drawparallels(numpy.arange(-90, 90, 30))
#
#    xh, yh = pmap(f[:,0], f[:,1])
#    l1 = pmap.scatter(xh, yh, 15, marker='o', c=f[:,4],
#                  edgecolor='none', vmin=-80, vmax=80)
#    pmap.colorbar(l1,location='bottom',pad="5%")
#    plt.title('Subducting Plate Velocity (km/Myr)')
#    plt.show()
#
#    pmap = Basemap(projection='hammer', lat_0=0, lon_0=180,
#               resolution='l')
#    pmap.drawmapboundary(fill_color='white')
#    pmap.fillcontinents(color='grey', lake_color='white', zorder=0)
#    pmap.drawmeridians(numpy.arange(0, 360, 30))
#    pmap.drawparallels(numpy.arange(-90, 90, 30))
#
#    xh, yh = pmap(f[:,0], f[:,1])
#    l1 = pmap.scatter(xh, yh, 15, marker='o', c=f[:,5],
#                  edgecolor='none', vmin=-50, vmax=50)
#    pmap.colorbar(l1,location='bottom',pad="5%")
#    plt.title('OP Velocity (km/Myr)')
#    plt.show()
#
#    pmap = Basemap(projection='hammer', lat_0=0, lon_0=180,
#               resolution='l')
#    pmap.drawmapboundary(fill_color='white')
#    pmap.fillcontinents(color='grey', lake_color='white', zorder=0)
#    pmap.drawmeridians(numpy.arange(0, 360, 30))
#    pmap.drawparallels(numpy.arange(-90, 90, 30))
#
#    xh, yh = pmap(f[:,0], f[:,1])
#    l1 = pmap.scatter(xh, yh, 15, marker='o', c=f[:,6],
#                  edgecolor='none', vmin=-50, vmax=50)
#    pmap.colorbar(l1,location='bottom',pad="5%")
#    plt.title('Trench Velocity (km/Myr)')
#    plt.show()
#
#    pmap = Basemap(projection='hammer', lat_0=0, lon_0=180,
#               resolution='l')
#    pmap.drawmapboundary(fill_color='white')
#    pmap.fillcontinents(color='grey', lake_color='white', zorder=0)
#    pmap.drawmeridians(numpy.arange(0, 360, 30))
#    pmap.drawparallels(numpy.arange(-90, 90, 30))
#
#    xh, yh = pmap(f[:,0], f[:,1])
#    l1 = pmap.scatter(xh, yh, 15, marker='o', c=f[:,7],
#                  edgecolor='none')
#    pmap.colorbar(l1,location='bottom',pad="5%")
#    plt.title('Subduction Obliquity (degrees)')
#    plt.show()
#
#    pmap = Basemap(projection='hammer', lat_0=0, lon_0=180,
#               resolution='l')
#    pmap.drawmapboundary(fill_color='white')
#    pmap.fillcontinents(color='grey', lake_color='white', zorder=0)
#    pmap.drawmeridians(numpy.arange(0, 360, 30))
#    pmap.drawparallels(numpy.arange(-90, 90, 30))
#
#    xh, yh = pmap(f[:,0], f[:,1])
#    l1 = pmap.scatter(xh, yh, 15, marker='o', c=f[:,8],
#                  edgecolor='none')
#    pmap.colorbar(l1,location='bottom',pad="5%")
#    plt.title('Subduction Polarity (degrees)')
#    plt.show()
#
#    pmap = Basemap(projection='hammer', lat_0=0, lon_0=180,
#               resolution='l')
#    pmap.drawmapboundary(fill_color='white')
#    pmap.fillcontinents(color='grey', lake_color='white', zorder=0)
#    pmap.drawmeridians(numpy.arange(0, 360, 30))
#    pmap.drawparallels(numpy.arange(-90, 90, 30))
#
#    xh, yh = pmap(f[:,0], f[:,1])
#    l1 = pmap.scatter(xh, yh, 15, marker='o', c=f[:,15],
#                  edgecolor='none', vmin=-50, vmax=1000)
#    pmap.colorbar(l1,location='bottom',pad="5%")
#    plt.title('Distance to Slab Edge (km)')
#    plt.show()
#
