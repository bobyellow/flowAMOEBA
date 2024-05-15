# Author and copyright information
__author__ = "Ran Tao"
__credits__ = "Copyright (c) 2018-01 Ran Tao"
__license__ = "New BSD License"
__version__ = "1.0.0"
__email__ = "rtao@usf.edu"

import time as tm
import numpy
import clusterpy
import operator
import csv
from componentsAlg import calculateGetisG, quickSort2, neighborSort, randomOD
from contiguity import weightsFromFlows  # obtain spatial weights of flow data

__all__ = ['execAMOEBA']


def execAMOEBAFLOW(AREAS1, AREAS2, FlowValue, significance, seed_threshold):
    """AMOEBA A Multidirectional Optimum Ecotope-Based Algorithm

    AMOEBA, devised by [Alstadt_Getis2006]_, embeds a local spatial
    autocorrelation statistic in an iterative procedure in order to
    identify spatial clusters (ecotopes) of related spatial units.


    This algorithm starts with an initial area to which neighboring areas are
    iteratively attached until the addition of any neighboring area fails to
    increase the magnitude of the local Getis statistic of [Getis_Ord1992]_
    and [Ord_Getis1995]_. The resulting region is considered an ecotope. This
    procedure is executed for all areas, and final ecotopes are defined after
    resolving overlaps and asserting non-randomness.


    The algorithm implemented here is a new version of AMOEBA proposed by
    [Duque_Alstadt_Velasquez_Franco_Betancourt2010]_ that significantly
    reduces its computational complexity without losing optimality. ::

        layer.cluster('amoeba',vars,<wType>,<significance>)

    :keyword vars: Area attribute(s)
    :type vars: sist
    :keyword wType: Type of first-order contiguity-based spatial matrix: 'rook' or 'quuun'. Default value wType = 'rook'.
    :type wType: string
    :keyword significance: Level of statistical significance. Default value
    :type significance: float (significance=0.01)

    IMPORTANT NOTE:

    Since AMEOBA is a non-exhausive algorithm, clusterPy does not provide the
    dissolve option. to obtain the solution vector you will need to export the
    layer with the command "Layer.exportArcData". The exported shape file will
    have an additional variable with the solution vector, where areas with
    positive values belongs to high value clusters; areas with negative values
    belongs to low value clusters; and areas with value zero are those outside
    the clusters.
    """
    start = tm.time()
    print "Running FlowAMOEBA by Ran Tao, derived from computationally efficient AMOEBA (Duque et al., 2010)"

    # Initialize the number of clusters as 0
    NumberOfClusters = 0
    print 'significance: ' + str(significance)
    areas1 = AREAS1  # list of O area
    areas2 = AREAS2  # list of D area
    flowvalue = FlowValue  ##OD pairs with non-zero value
    y1 = areas1.Y  # Read OD features
    y2 = areas2.Y
    Olength = len(y1)
    Dlength = len(y2)
    Flowlength = Olength * Dlength

    # Read flow data
    y = flowvalue
    yKeys = y.keys()
    # print len(yKeys)
    for i in range(len(yKeys)):  # remove zero-value flow
        if yKeys[i][0] == 0 or yKeys[i][1] == 0:
            del y[yKeys[i]]
    # Get the flow dictionary keys, i.e. (O,D) tuple
    yKeys = y.keys()
    # print len(yKeys)

    areaKeys = y.keys()  ##y is a dictionary of input polygons: key+values

    # Obtain flow neighbors based on O&D contiguity, by default Rook's case
    # Wflow is a dictionary, for each item (flow), the key is its (O,D) tuple, the value includes all the (O,D) tuples of its neighbor flows

    Wflow = weightsFromFlows(areas1, areas2,y,1)
    # the last parameter decides the level of flow neighbors.
    # level 1 has O (D) the same, while D (O) is neighbor; level 2 has both O and D as neighbor; level 12 = level 1 + level 2


    # wValues = Wflow.values()
    # print Wflow[(9,183)]
    w = Wflow  # flow weights: binary

    time1 = tm.time() - start
    print 'finish calculating Wflow, running time: ' + str(time1) + ' seconds'

    # Calculate mean, sum, std for calculating G* later
    dataLength = len(y)
    dataSum = numpy.sum(numpy.double(y.values()))
    dataMean = numpy.mean(numpy.double(y.values()))
    dataStd = numpy.std(numpy.double(y.values()))
    print 'dataLength: ' + str(dataLength)
    print 'dataSum: ' + str(dataSum)
    print 'dataMean:' + str(dataMean)
    print 'dataStd: ' + str(dataStd)

    # Define cluster and its G* as dictionary
    generatedClusters = dict()
    clusterGValues = dict()
    clusterGValuesAbs = dict()
    #
    # print "Start iterative process"

    # Iterate every flow as seed
    for s in areaKeys:
        discNeighbor = []
        if y[s] < seed_threshold:
            continue
        if s in w:  # w #w = Wflow, dictionary of flow's neighbors
            neighbors = w[s]  # Get the neighbors of flow s
        itAreaList = [s]

        # calculate local G* for each flow, then calculate local G for itself and its neighbors
        currentG = calculateGetisG([s], dataMean, dataStd, y, dataLength)
        previousG = currentG - 1
        sortedNeighbors = []
        while currentG != previousG:
            # Sort the neighbors by its flow value
            NeighborsDic = {x: y[x] for x in neighbors}  # get neighbor as dictionary
            sortedNeighbors = []
            sortedNeighborsDic = sorted(NeighborsDic.items(), key=operator.itemgetter(1))
            for nei in sortedNeighborsDic:
                sortedNeighbors.append(nei[0])
            neighbors = []
            previousG = currentG
            AuxItAreaList = []
            AuxDiscNeighbor = []
            AreasBase = itAreaList
            NoSorted = len(sortedNeighbors)  # Number of sorted neighbors

            if currentG <= 0:  # Cold spot
                for a in range(
                        NoSorted):  # Try to include seed flow's neighbor one by one, from the largest to the smallest
                    newG = calculateGetisG(itAreaList + sortedNeighbors[0: a + 1], dataMean, dataStd, y, dataLength)
                    if newG < currentG:  ##even colder...include this neighbor into the ecotope
                        currentG = newG  # update the G* value after including a new flow
                        AuxItAreaList = sortedNeighbors[0: a + 1]  # add the neighbor to the ecotope
                        AuxDiscNeighbor = sortedNeighbors[
                                          a + 1: NoSorted]  # remove the neighbor from the sortedNeighbor list
            else:  # Hot spot
                for a in range(NoSorted):
                    newG = calculateGetisG(itAreaList + sortedNeighbors[-a - 1: NoSorted], dataMean, dataStd, y,
                                           dataLength)
                    if newG > currentG:  ##even hotter....include this neighbor into the ecotope
                        currentG = newG
                        AuxItAreaList = sortedNeighbors[-a - 1: NoSorted]
                        AuxDiscNeighbor = sortedNeighbors[0: -a - 1]

            discNeighbor = discNeighbor + AuxDiscNeighbor
            itAreaList = itAreaList + AuxItAreaList
            for x in AuxItAreaList:
                neighbors = neighbors + list(
                    set(w[x]) - set(sortedNeighbors) - set(itAreaList) - set(discNeighbor) - set(neighbors))

        # Save the cluster (a list of (O,D) keys and G* value) grown from this seed flow
        generatedClusters[s] = itAreaList
        clusterGValues[s] = currentG
        clusterGValuesAbs[s] = numpy.abs(currentG)
        # continue to the next seed flow

    neighborSorted = []
    neighborSortedDic = sorted(clusterGValuesAbs.items(), key=operator.itemgetter(1))
    for nei in neighborSortedDic:
        neighborSorted.append(nei[0])

    prioritaryClusters = reversed(neighborSorted)

    tnclusterCounter = 1
    output = {}
    clusterMap = {}
    mapCounter = 0
    areaRange = range(dataLength)
    # areaRange = range(dataNonZeroLength)
    randomKeyList = []
    clusterCounter = 0
    #####################################################################

    time = tm.time() - start
    print 'finish calculating FlowAMOEBA for observations, running time: ' + str(time) + ' seconds'
    print "Start testing clusters significance"
    ##generatedClusters only includes clustered (pos/neg) objects and their clustering neighbor
    ##output stores the AMOEBA Value
    for i in areaKeys:
        output[i] = 0
        clusterMap[i] = mapCounter
        mapCounter = mapCounter + 1
    negClusterCounter = -1
    posClusterCounter = 1

    # Original Permutation Process#############################################
    # the original permutation creates a random list of Keys without changing flow value, basically it randomly switches flow (O,D) key

    for i in range(1000):
        randomKeyListInstance = []
        randomList = numpy.random.permutation(areaRange)  ##reshuffle Flow ID, keep total number the same
        for j in randomList:
            randomKeyListInstance.append(areaKeys[j])
        randomKeyList.append(randomKeyListInstance)  ##create a randomly permuted ID(key) list

    for x in prioritaryClusters:  # Permutation and overlapping tests
        validCluster = 1
        for h in generatedClusters[x]:  # Validates if the cluster overlaps with a previously accepted one
            if output[h] != 0:
                validCluster = 0
                break
        if validCluster == 1:
            betterClusters = 0
            for j in range(1000):  # Monte Carlo permutation test
                permKey = []
                for e in generatedClusters[x]:
                    t = clusterMap[e]
                    permKey.append((randomKeyList[j])[t])
                ## permKey is a randomly generated list of flow key, list length = cluster length
                randomG = calculateGetisG(permKey, dataMean, dataStd, y, dataLength)
                # print randomG
                if clusterGValues[x] >= 0:
                    if clusterGValues[x] < randomG:
                        betterClusters = betterClusters + 1
                else:
                    if clusterGValues[x] > randomG:
                        betterClusters = betterClusters + 1
            pValue = (betterClusters) / 1000.00
            if pValue <= significance:
                NumberOfClusters = NumberOfClusters + 1
                clustId = 0
                if numpy.double(clusterGValues[x]) > 0:
                    clustId = posClusterCounter
                    posClusterCounter = posClusterCounter + 1
                else:
                    clustId = negClusterCounter
                    negClusterCounter = negClusterCounter - 1
                for h in generatedClusters[x]:
                    if clusterGValues[x] < 0:
                        output[h] = clustId
                    else:
                        output[h] = clustId
                clusterCounter = clusterCounter + 1
    ############  END OF SIMULATION  ###############################

    print 'number of positive clusters: ' + str(posClusterCounter - 1)
    print 'number of negative clusters: ' + str(abs(negClusterCounter) - 1)
    print 'number of total clusters: ' + str(clusterCounter)

    time2 = tm.time() - start
    time2nd = time2 - time1
    print '2nd half running time: ' + str(time2nd) + ' seconds'

    time0 = tm.time() - start
    print 'total running time: ' + str(time0) + ' seconds'

    outputStr = 'O, D, AMOEBA, FlowValue'
    for i in output:
        outputStr = outputStr + '\n' + str(i) + ', ' + str(output[i]) + ', ' + str(y[i])

    outputStr1 = outputStr.replace('(', '')
    outputStr2 = outputStr1.replace(')', '')

    return outputStr2
