# -*- coding: utf-8 -*-
"""
MFIX [Multiphase Flow with Interphase eXchanges] is a general-purpose
computer code developed at the National Energy Technology Laboratory
[NETL] for describing the hydrodynamics, heat transfer and chemical
reactions in fluid-solid systems.

Please visit: https://mfix.netl.doe.gov/

This python file contains code for parsing, fixing, and writing stl files

Last update: 1/31/2014

@author: Justin Weber
"""
# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

import logging
import re
import struct
import math
import os

mylogger = logging.getLogger(__name__)
try:
    from mfixgui.guilib import config
    CONFIG = config.guiConfig()
    mylogger.setLevel(CONFIG.get('debugging', 'level').upper())
except:
    pass

class StlParser(object):
    def __init__(self, stlFile = None):
        """
        Parse an ASCII or Binary STL file.

        Usage:
        stlparser = StlParser(stlFile = './path/to/stl/file/fname.stl')
        """

        # Variables
        self.stlFileType = None
        self.solids = {}
        self.stats = {}
        self.changesMade = False
        self.solidNumDict={}
        self.avaliableFacets = []

        if stlFile:
            self.readStlFile(stlFile = stlFile)

    # --- stl read/write functions ---
    def saveStlFile(self, fname, ftype = 'ASCII'):
        """
        Save the stl data to either an ASCII or Binary stl file.

        Usage:
        saveStlFile('filename', ftype = 'ASCII')
        """
        if fname:
            # place unused facets in another solid

            if self.solidNumDict.keys():
                num=min(set(range(1,max(self.solidNumDict.keys())+2))-set(self.solidNumDict.keys()))
            else:
                num=1

            name = 'solid_{}'.format(num)
            if name not in self.solids:
                self.solids[name] = self.avaliableFacets
                self.solidNumDict[num] = name
            else:
                self.solids[name+'_2'] = self.avaliableFacets
                self.solidNumDict[num] = name+'_2'

            self.solids.pop('avaliable')
            self.solidNumDict

            if ftype.lower() == 'ascii':
                self.writeAsciiStlFile(fname)
            elif ftype.lower() == 'binary' or ftype.lower() == 'b':
                self.writeBinaryStlFile(fname)
            else:
                mylogger.warning('Unknown stl file type: %s', ftype)
                pass
        else:
            mylogger.debug('No stl filename specified.')
            pass

    def writeAsciiStlFile(self, fname):
        """
        Write ASCII stl file
        """
        mylogger.debug('Writing ASCII stl file: %s', fname)

        with open(fname, 'w') as stlFile:
            for solid in self.solidNumDict.keys():
                if len(self.solids[self.solidNumDict[solid]])>0:
                    stlFile.write('solid {0:04d}\n'.format(solid))

                    for facet in self.solids[self.solidNumDict[solid]]:
                        stlFile.write('facet normal {0:e} {1:e} {2:e}\n'.format(facet[0][0], facet[0][1], facet[0][2]))
                        stlFile.write('    outer loop\n')
                        for point in facet[1]:
                            stlFile.write('        vertex {0:e} {1:e} {2:e}\n'.format(point[0], point[1], point[2]))
                        stlFile.write('    endloop\n')
                        stlFile.write('endfacet\n')
                    stlFile.write('endsolid {0:04d}\n'.format(solid))

    def writeBinaryStlFile(self, fname):
        """
        Write Binary stl file
        """
        mylogger.debug('Writing Binary stl file: %s', fname)

        with open(fname, 'wb') as stlFile:
            stlFile.write(struct.pack(b'80sI', b'MFIX GUI STL Writer', self.stats['facets']))
            for solid in self.solidNumDict.keys():
                if len(self.solids[self.solidNumDict[solid]])>0:

                    for facet in self.solids[self.solidNumDict[solid]]:
                        facetData = [facet[0][0], facet[0][1], facet[0][2]]
                        for point in facet[1]:
                            facetData+=[point[0], point[1], point[2]]
                        facetData.append(0.0)

                        stlFile.write(struct.pack(b'12fH', *facetData))

    def readStlFile(self, stlFile = None):
        """
        Parse an ASCII or Binary STL file.
        """
        if not stlFile:
            return

        self.stlFile = stlFile

        mylogger.debug('Parsing STL File: %s', str(self.stlFile))

        if self.isStlAscii():
            self.stlFileType = 'ASCII'
            self.stats['type'] = 'ASCII'
            self.parseAsciiFile()
        else:
            self.stlFileType = 'Binary'
            self.stats['type'] = 'Binary'
            self.parseBinaryFile()

        self.genFacetLists()
        self.calculateStats()

    def parseAsciiFile(self):
        """
        Read a ASCII formated STL file.
        """

        # Regex
        solidRegex = re.compile("solid(.+?)endsolid", re.IGNORECASE|re.DOTALL|re.MULTILINE)
        facetRegex = re.compile("facet(.+?)endfacet", re.IGNORECASE|re.DOTALL|re.MULTILINE)
        loopRegex = re.compile("outer loop(.+?)endloop", re.IGNORECASE|re.DOTALL|re.MULTILINE)
        floatRegex = re.compile('([+-]?\d+(?:\.\d+)(?:[eE][+-]?\d*)?)')
        normalRegex = re.compile("normal *([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?) *([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?) *([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", re.IGNORECASE|re.DOTALL|re.MULTILINE)
        vertexRegex = re.compile("vertex *([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?) *([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?) *([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", re.IGNORECASE|re.DOTALL|re.MULTILINE)

        facetTotal = 0

        with open(self.stlFile, 'rU') as stlFile:
            stlFileLines = stlFile.read()

        # Look for solids
        for solidCount, solid in enumerate(solidRegex.findall(stlFileLines)):
            # Look for solid name
            solidName = solid.split('\n')[0].strip()
            if solidName:
                num = re.findall('(\d+)', solidName)
                if num:
                    num = int(num[-1])
                    if num==0 or num in self.solidNumDict:
                        num=min(set(range(1,max(list(self.solidNumDict.keys())+[1])+2))-set(self.solidNumDict.keys()))
                else:
                    num=min(set(range(1,max(list(self.solidNumDict.keys())+[1])+2))-set(self.solidNumDict.keys()))

                if not solidName in self.solids:
                    self.solids[solidName]=[]
                    self.solidNumDict[int(num)]=solidName
                else:
                    print('duplicate solid')

            # Look for facets
            for facetCount, facet in enumerate(facetRegex.findall(solid)):

                # parse facet
                vertexList = []
                normal = []
                loop = False
                loopCount = 0
                vertexCount = 0
                parseNormal = False
                parseVertex = False
                for line in facet.split('\n'):

                    if 'normal' in line:
                        parseNormal = True
                    if parseNormal:
                        normal += [float(x) for x in floatRegex.findall(line)]

                    if 'outer loop' in line:
                        parseNormal = False
                        loop = True
                        loopCount+=1
                    elif 'endloop' in line:
                        loop = False

                    if 'vertex' in line:
                        parseVertex = True
                    if loop and parseVertex:
                        vertexList += [float(x) for x in floatRegex.findall(line)]

                # check and split vertexes
                vertexCount = int(len(vertexList)/3.0)
                vertexList = [tuple(vertexList[i:i+3]) for i in range(0, len(vertexList), 3)]

                # Do some checking
                if loopCount > 1:
#                    mylogger.warning('Multiple facets loops: %s', str(facetCount+1))
                    pass

                if vertexCount > 3:
#                    mylogger.warning('Found more than three vertexes: %s', str(facetCount+1))
                    pass

                self.solids[solidName].append((tuple(normal), tuple(vertexList)))

            facetTotal+=(facetCount+1)
            self.stats[solidName] = {'facets':facetCount+1}

#        mylogger.debug('Parsed %s solid(s)', str(solidCount+1))
        self.stats['solids'] = solidCount+1
#        mylogger.debug('Parsed %s facet(s)', str(facetTotal))
        self.stats['totalfacetsRead'] = facetTotal


    def parseBinaryFile(self):
        """
        Read a Binary STL file.
        """
        headerLength = 80
        facetCountLength = 4
        facetLength = 12*4 + 2 #twelve 32 bit floats and one 2 bit unsigned int

        self.solids = {'solid1':[]}

        with open(self.stlFile, 'rb') as stlFile:
            stlFile.seek(headerLength)
            facetTotal = struct.unpack('<I',stlFile.read(facetCountLength))[0]
            for facet in range(facetTotal):
                raw = stlFile.read(facetLength)
#                try:
                facetData = struct.unpack('<12fH', raw)
#                except:
#                    mylogger.debug('bad facet data, len: {}'.format(len(raw)))
#                    continue

                vertexList = []

                for i in range(3,12,3):
                    vertexList.append((facetData[i],
                                       facetData[i+1],
                                       facetData[i+2]))

                normal = (facetData[0], facetData[1], facetData[2])

                self.solids['solid1'].append((normal, tuple(vertexList)))

        mylogger.debug('Parsed %s solid(s)', str(1))
        self.stats['solids'] = 1
        mylogger.debug('Parsed %s facet(s)', str(facetTotal))
        self.stats['totalfacetsRead'] = facetTotal


    def isStlAscii(self):
        """
        Check if stl file is ASCII

        A proper ASCII stl file should begin with "solid" and end with "endsolid"
        """

        with open(self.stlFile,'rb') as stlFile:
            # Check first line if ASCII, should start with "solid"
            try:
                solid = stlFile.readline().strip().lower().startswith(b'solid')
            except:
                solid = False


            # Double check by reading triangle count if binary
            stlFile.seek(0)
            stlFile.seek(80)
            facetTotal = struct.unpack(b'<I', stlFile.read(4))[0]


            # Check Last Line
            if solid:
                stlFile.seek(-2, 2)
                while stlFile.read(1) != b"\n":
                    stlFile.seek(-2, 1)
                endsolid = stlFile.readline().strip().lower().startswith(b'endsolid')


        # Check if facetTotal is correct number by comparing to the filesize
        if facetTotal * 50 + 84 == os.path.getsize(self.stlFile):
            solid = False

        if not solid:
            mylogger.debug('Binary STL File Detected')
            return False
        else:
            mylogger.debug('ASCII STL File Detected')
            if endsolid:
                return True
            else:
                mylogger.warning('No "endsolid" detected: poorly formated, or missing data')
                return True

    def genFacetLists(self):

        self.avaliableFacets=[]
        removelist = []
        for solid in self.solids:
            self.avaliableFacets+=self.solids[solid]

            for key, value in list(self.solidNumDict.items()):
                if value == solid:
                    self.solidNumDict.pop(key)
            removelist.append(solid)

        for solid in removelist:
            self.solids.pop(solid)
        self.solids={'avaliable':self.avaliableFacets}


    # --- STL error checking and cleanup ----
    def cleanStlFile(self):
        """
        Try to clean the stl file of duplicate faces and degenerate faces.
        """
        self.checkForDegenerateFacets()
        self.checkForDuplicateFacets()

    def checkForDegenerateFacets(self):
        """
        Look for all facets that are not triangles.

        A traingle is defined with three unique points.
        """
        for solid in self.solids:
            removeList = []
            for facet in self.solids[solid]:
                if len(set(facet[1]))!=3:
                    removeList.append(facet)
                    self.changesMade = True

            mylogger.debug('Removed %s degenerate facets', str(len(removeList)))
            self.stats[solid]['degeneratefacets'] = len(removeList)
            self.solids[solid][:] = [i for i in self.solids[solid] if not i in removeList]

    def checkForDuplicateFacets(self):
        """
        look for identical copies of facets: (normal,(pt1, pt2, pt3))
        """
        for solid in self.solids:
            uniqueFacets = set(self.solids[solid])
            if len(self.solids[solid]) != len(uniqueFacets):
                mylogger.debug('Removed %s duplicate facets', str(len(self.solids[solid]) - len(uniqueFacets)))
                self.stats[solid]['duplicatefacets'] = len(self.solids[solid]) - len(uniqueFacets)
                self.solids[solid] = list(uniqueFacets)
                self.changesMade = True
            else:
                self.stats[solid]['duplicatefacets'] = 0

    def calculateNewNormalVectors(self):
        """
        calculate new normal vectors based on the right hand rule.
        """
        for solid in self.solids:
            for i, facet in enumerate(self.solids[solid]):
                unitVect = self.calcUnitVect(facet[1])

                self.solids[solid][i] = (unitVect, facet[1])

    # --- clalculations ---
    def calcUnitVect(self, points):
        u = self._ptSub(points[1], points[0])
        v = self._ptSub(points[2], points[0])

        uVect = (u[1]*v[2] - u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])

        return uVect

    def _ptSub(self, p1, p2):
        sol = []
        for i in range(len(p1)):
            sol.append(p1[i]-p2[i])
        return tuple(sol)


    def dotproduct(self, v1, v2):
        return sum((a*b) for a, b in zip(v1, v2))

    def length(self, v):
        return math.sqrt(self.dotproduct(v, v))

    def angle(self, v1, v2):
        return math.acos(self.dotproduct(v1, v2) / (self.length(v1) * self.length(v2)))*180.0/math.pi

    def pointinbox(self, point, x, y, z):
        if point[0]>=x[0] and point[0]<=x[1] and \
                point[1]>=y[0] and point[1]<=y[1] and \
                point[2]>=z[0] and point[2]<=z[1]:
            return True
        else:
            return False

    def pointinellipsoid(self, point, x, y, z):

        center = [sum(axis)/2.0 for axis in [x, y, z]]
        axes = [max(axis)-c for axis, c in zip([x, y, z], center)]

        if any(axis==0 for axis in axes):
            return False
        if sum([((p-offset)/axis)**2 for p, axis, offset in zip(point, axes, center)])<=1:
            return True
        else:
            return False

    def facetinbox(self, facet, x, y, z, greedy):
        boolList=[]
        for point in facet[1]:
            boolList.append(self.pointinbox(point, x, y, z))

        if greedy:
            return any(boolList)
        elif not greedy:
            return all(boolList)
        else:
            return False

    def facetinellipsoid(self, facet, x, y, z, greedy):
        boolList=[]
        for point in facet[1]:
            boolList.append(self.pointinellipsoid(point, x, y, z))

        if greedy:
            return any(boolList)
        elif not greedy:
            return all(boolList)
        else:
            return False

    def invertNormal(self, facet):
        points = [facet[1][-1], facet[1][1], facet[1][0]]

        uvect = self.calcUnitVect(points)

        return [uvect, points]


    # --- stats and report data ---
    def calculateStats(self):
        """
        loop through all facets, finding extents, aspect ratios
        """

        totalaspectRatioList = []
        totalfacetSegmentList = []

        for solid, facets in list(self.solids.items()):

            if len(facets)==0:
                self.solids.pop(solid)
                if solid in self.stats:
                    self.stats.pop(solid)
                continue

            minmaxdict = {'xMin':float('+inf'),'xMax':float('-inf'),
                          'yMin':float('+inf'),'yMax':float('-inf'),
                          'zMin':float('+inf'),'zMax':float('-inf')}

            if solid in self.stats:
                self.stats[solid].update(minmaxdict)
            else:
                self.stats[solid] = minmaxdict

            aspectRatioList = []
            facetSegmentList = []

            for facet in facets:
                for point in facet[1]:
                    # Look for extents
                    self.stats[solid]['xMin'] = min(self.stats[solid]['xMin'], point[0])
                    self.stats[solid]['xMax'] = max(self.stats[solid]['xMax'], point[0])
                    self.stats[solid]['yMin'] = min(self.stats[solid]['yMin'], point[1])
                    self.stats[solid]['yMax'] = max(self.stats[solid]['yMax'], point[1])
                    self.stats[solid]['zMin'] = min(self.stats[solid]['zMin'], point[2])
                    self.stats[solid]['zMax'] = max(self.stats[solid]['zMax'], point[2])

                # Calculate facet line lengths
                a = ((facet[1][0][0]-facet[1][1][0])**2 + (facet[1][0][1]-facet[1][1][1])**2 + (facet[1][0][2]-facet[1][1][2])**2)**(0.5)
                b = ((facet[1][1][0]-facet[1][2][0])**2 + (facet[1][1][1]-facet[1][2][1])**2 + (facet[1][1][2]-facet[1][2][2])**2)**(0.5)
                c = ((facet[1][2][0]-facet[1][0][0])**2 + (facet[1][2][1]-facet[1][0][1])**2 + (facet[1][2][2]-facet[1][0][2])**2)**(0.5)

                facetSegmentList+=[a, b, c]

                # Calculate facet aspect ratio: circumradius/(2*inradius)
                s = (a + b + c) / 2.0
                try:
                    aspectRatio = (a * b * c) / (8.0 * (s - a) * (s - b) * (s - c))
                    aspectRatioList.append(aspectRatio)
                except ZeroDivisionError:
                    pass

            totalaspectRatioList+=aspectRatioList
            totalfacetSegmentList+=facetSegmentList

            self.stats[solid]['facets'] = len(facets)
            self.stats[solid]['width'] = self.stats[solid]['xMax']-self.stats[solid]['xMin']
            self.stats[solid]['height'] = self.stats[solid]['yMax']-self.stats[solid]['yMin']
            self.stats[solid]['depth'] = self.stats[solid]['zMax']-self.stats[solid]['zMin']
            self.stats[solid]['aspectRatio'] = {'min':min(aspectRatioList),'max':max(aspectRatioList),'average':sum(aspectRatioList)/len(aspectRatioList)}
            self.stats[solid]['facetSegment'] = {'min':min(facetSegmentList),'max':max(facetSegmentList),'average':sum(facetSegmentList)/len(facetSegmentList)}

        # Update Report Data
        self.stats['aspectRatio'] = {'min':min(totalaspectRatioList),'max':max(totalaspectRatioList),'average':sum(totalaspectRatioList)/len(totalaspectRatioList)}
        self.stats['facetSegment'] = {'min':min(totalfacetSegmentList),'max':max(totalfacetSegmentList),'average':sum(totalfacetSegmentList)/len(totalfacetSegmentList)}
        self.stats['xMin'] = min([self.stats[solid]['xMin'] for solid in self.stats if isinstance(self.stats[solid], dict) and 'xMin' in self.stats[solid]])
        self.stats['xMax'] = max([self.stats[solid]['xMax'] for solid in self.stats if isinstance(self.stats[solid], dict) and 'xMax' in self.stats[solid]])
        self.stats['yMin'] = min([self.stats[solid]['yMin'] for solid in self.stats if isinstance(self.stats[solid], dict) and 'yMin' in self.stats[solid]])
        self.stats['yMax'] = max([self.stats[solid]['yMax'] for solid in self.stats if isinstance(self.stats[solid], dict) and 'yMax' in self.stats[solid]])
        self.stats['zMin'] = min([self.stats[solid]['zMin'] for solid in self.stats if isinstance(self.stats[solid], dict) and 'zMin' in self.stats[solid]])
        self.stats['zMax'] = max([self.stats[solid]['zMax'] for solid in self.stats if isinstance(self.stats[solid], dict) and 'zMax' in self.stats[solid]])
        self.stats['facets'] = sum([self.stats[solid]['facets'] for solid in self.stats if isinstance(self.stats[solid], dict) and 'facets' in self.stats[solid]])
        self.stats['width'] = self.stats['xMax']-self.stats['xMin']
        self.stats['height'] = self.stats['yMax']-self.stats['yMin']
        self.stats['depth'] = self.stats['zMax']-self.stats['zMin']
        self.stats['center'] = [(self.stats['xMin']+self.stats['xMax'])/2.0,
                                (self.stats['yMin']+self.stats['yMax'])/2.0,
                                (self.stats['zMin']+self.stats['zMax'])/2.0,
                                ]

        return self.stats

    # --- stl modifcation functions ---
    def scale(self, axes, factor):
        """
        Scale the stl file in the provided axes by the factor about the center of the stl file.

        Usage:
        scale('xyz', 4)
        """
        if axes.lower() == 'all':
            axes = 'xyz'

        center  = self.stats['center']

        for solid in self.solids:
            for i, facet in enumerate(self.solids[solid]):
                newFacet = [facet[0],list(facet[1])]
                for j, point in enumerate(facet[1]):
                    newPoint = list(point)
                    if 'x' in axes.lower():
                        newPoint[0] = (point[0]-center[0])*factor+center[0]
                    if 'y' in axes.lower():
                        newPoint[1] = (point[1]-center[1])*factor+center[1]
                    if 'z' in axes.lower():
                        newPoint[2] = (point[2]-center[2])*factor+center[2]
                    newFacet[1][j] = tuple(newPoint)
                self.solids[solid][i] = (newFacet[0], tuple(newFacet[1]))
        self.calculateNewNormalVectors()
        self.calculateStats()

    def translate(self, axes, factor):
        """
        Transform the stl file in the provided axes by the factor

        Usage:
        transform('xyz', 4)
        """
        if axes.lower() == 'all':
            axes = 'xyz'

        for solid in self.solids:
            for i, facet in enumerate(self.solids[solid]):
                newFacet = [facet[0],list(facet[1])]
                for j, point in enumerate(facet[1]):
                    newPoint = list(point)
                    if 'x' in axes.lower():
                        newPoint[0] = point[0]+factor
                    if 'y' in axes.lower():
                        newPoint[1] = point[1]+factor
                    if 'z' in axes.lower():
                        newPoint[2] = point[2]+factor
                    newFacet[1][j] = tuple(newPoint)
                self.solids[solid][i] = (newFacet[0], tuple(newFacet[1]))
        self.calculateNewNormalVectors()
        self.calculateStats()

    def rotate(self, axis, degrees):
        """
        rotate the stl file around the provided axis by the specified degrees

        Usage:
        rotate('x', 45)
        """
        theta = degrees/180.0*math.pi

        for solid in self.solids:
            for i, facet in enumerate(self.solids[solid]):
                newFacet = [facet[0],list(facet[1])]
                for j, point in enumerate(facet[1]):
                    newPoint = list(point)
                    if 'x' == axis.lower():
                        newPoint[2] = point[2]*math.cos(theta)-point[1]*math.sin(theta)
                        newPoint[1] = point[2]*math.sin(theta)+point[1]*math.cos(theta)
                    elif 'y' == axis.lower():
                        newPoint[0] = point[0]*math.cos(theta)-point[2]*math.sin(theta)
                        newPoint[2] = point[0]*math.sin(theta)+point[2]*math.cos(theta)
                    elif 'z' == axis.lower():
                        newPoint[0] = point[0]*math.cos(theta)-point[1]*math.sin(theta)
                        newPoint[1] = point[0]*math.sin(theta)+point[1]*math.cos(theta)
                    else:
                        mylogger.warning('could not determine rotation axis from: {}'.format(axis))
                        pass
                    newFacet[1][j] = tuple(newPoint)
                self.solids[solid][i] = (newFacet[0], tuple(newFacet[1]))

        self.calculateNewNormalVectors()
        self.calculateStats()


    # --- selections ---
    def newRegion(self, x=[0,0], y=[0,0], z=[0,0], vect=None, devAngle=0.0,
                  greedy=False, num=None, vol='box'):
        """
        Find facets that are inside a box or ellipsoid.
        """

        newSolid = []
        removelist = []
        for facet in self.avaliableFacets:
            if self.facetinbox(facet, x,y,z, greedy) and vol=='box' or\
                    self.facetinellipsoid(facet, x, y, z, greedy) and vol=='ellipsoid':
                if vect:
                    if abs(self.angle(facet[0], vect))<devAngle:
                        newSolid.append(facet)
                        removelist.append(facet)
                else:
                    newSolid.append(facet)
                    removelist.append(facet)

        for facet in removelist:
            self.avaliableFacets.remove(facet)

        if not num or num in self.solidNumDict:
            num=min(set(range(1,max(self.solidNumDict.keys())+2))-set(self.solidNumDict.keys()))

        name = 'solid_{}'.format(num)
        if name not in self.solids:
            self.solids[name] = newSolid
            self.solidNumDict[num] = name
        else:
            self.solids[name+'_2'] = newSolid
            self.solidNumDict[num] = name+'_2'

    def deleteRegion(self, x=[0,0], y=[0,0], z=[0,0], vect=None, devAngle=0.0,
                     greedy=False):

        for solid in self.solids:
            for facet in self.solids[solid][:]:
                if self.facetinbox(facet, x,y,z, greedy):
                    if vect:
                        if abs(self.angle(facet[0], vect))<devAngle:
                            self.solids[solid].remove(facet)
                    else:
                        self.solids[solid].remove(facet)

        removelist = []
        for facet in self.avaliableFacets:
            if self.facetinbox(facet, x,y,z, greedy):
                if vect:
                    if abs(self.angle(facet[0], vect))<devAngle:
                        removelist.append(facet)
                else:
                    removelist.append(facet)

        for facet in removelist:
            self.avaliableFacets.remove(facet)

    def invertNormalRegion(self, x=[0,0], y=[0,0], z=[0,0], vect=None, devAngle=0.0,
                           greedy=False):

        for solid in self.solids:
            for i, facet in enumerate(self.solids[solid]):
                if self.facetinbox(facet, x,y,z, greedy):
                    if vect:
                        if abs(self.angle(facet[0], vect))<devAngle:
                            self.solids[solid][i] = self.invertNormal(facet)
                    else:
                        self.solids[solid][i] = self.invertNormal(facet)

        for i, facet in enumerate(self.avaliableFacets):
            if self.facetinbox(facet, x,y,z, greedy):
                if vect:
                    if abs(self.angle(facet[0], vect))<devAngle:
                        self.avaliableFacets[i] = self.invertNormal(facet)
                else:
                    self.avaliableFacets[i] = self.invertNormal(facet)

def test_bunny(path):
    fname = 'stanfordbunny.stl'
    stlParser = StlParser(os.path.join(path, fname))

    for key, value in stlParser.stats.items():
        print('{} : {}'.format(key, value))

    #stlParser.deleteRegion([-20, 0], [-30,-10], [20, 40], greedy=True)
    #stlParser.invertNormalRegion([-35, -25], [-30,-20], [25, 35], greedy=True)
    stlParser.newRegion([-20, 0], [-30,-10], [20, 40], num=1, greedy=False, vol='box') # -10, -20, 30, d=20
    stlParser.newRegion([5, 25], [-30,-10], [20, 40], num=2, greedy=True, vol='box') # 15, -20, 30, d=20
    stlParser.newRegion([-45, -25], [-30,-10], [20, 40], num=3, greedy=True, vol='box', vect=[0,0,-1.0], devAngle=70) # -35, -20, 30, d=20
    stlParser.newRegion([-5, 5], [5,15], [75, 85], num=4, greedy=False, vol='ellipsoid') # 0, 10, 80, r=5
    stlParser.newRegion([0, 20], [-25,-5], [40, 60], num=5, greedy=False, vol='ellipsoid') # 10, -15, 50, r=10

    stlParser.saveStlFile(os.path.join(path, 'test.stl'), ftype='ascii')
    print('Done')

def test_sphere(path):
    fname = 'icosphere2.stl'
    stlParser = StlParser(os.path.join(path, fname))

    for key, value in stlParser.stats.items():
        print('{} : {}'.format(key, value))

    stlParser.newRegion([0, 2], [-1,1], [-1, 1], num=2, greedy=True, vol='ellipsoid')

    stlParser.saveStlFile(os.path.join(path, 'test.stl'), ftype='ascii')

    print('Done')

def test_flipnormal(path):
    fname = 'cylinder2.stl'

    stlParser = StlParser(os.path.join(path, fname))

    for key, value in stlParser.stats.items():
        print('{} : {}'.format(key, value))

    stlParser.invertNormalRegion([stlParser.stats['xMin'], stlParser.stats['xMax']],
                                 [stlParser.stats['yMin'], stlParser.stats['yMax']],
                                 [stlParser.stats['zMin'], stlParser.stats['zMax']]
                                 )

    stlParser.saveStlFile(os.path.join(path, fname), ftype='b')

    print('Done')

def test_scale(path):
    fname = 'icosphere2.stl'

    stlParser = StlParser(os.path.join(path, fname))

    for key in ['xMin', 'xMax']:
        print('{} : {}'.format(key, stlParser.stats[key]))

    stlParser.translate('xyz', 1)
    stlParser.scale('xyz', 4)

    for key in ['xMin', 'xMax']:
        print('{} : {}'.format(key, stlParser.stats[key]))

    stlParser.saveStlFile(os.path.join(path, 'test.stl'), ftype='a')


    print('Done')

def test_read():

    fname = 'icosphere2.stl'
    stlParser = StlParser(os.path.join(path, fname))

    for key, value in stlParser.stats.items():
        print('{} : {}'.format(key, value))

if __name__ == "__main__":

    path = os.path.abspath(os.path.join(os.path.abspath(__file__), os.pardir, os.pardir, ))
    path = os.path.join(path, 'data', 'stlfiles')

    #test_scale(path)
    #test_flipnormal(path)
    test_bunny(path)
#    test_sphere(path)
#    test_read()
