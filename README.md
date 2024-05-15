# flowAMOEBA
flowAMOEBA: a bottom-up flow clustering algorithm 
Ran Tao, Jean-Claude Thill. (2019). flowAMOEBA: Identifying Regions of Anomalous Spatial Interactions. _Geographical Analysis_ 51(1): 111–130.

Before running the code, there are several steps you need to prepare. First, import the origin and destination shapefile.
First, import the origin and destination shapefile.

AREAS1 = clusterpy.importArcData("yourpath/Origin_shapefile")
AREAS2 = clusterpy.importArcData("yourpath/Destination_shapefile")

Second, this is the execute code for FlowLISA
execAMOEBAFLOW(AREAS1, AREAS2, FlowValue, significance, seed_threshold)
    Parameters:
    - AREAS1: List of Origin areas.
    - AREAS2: List of Destination areas.
    - FlowValue: OD pairs with non-zero value.
    - Significance: Specifies the significance level.
    - Seed_threshold: Set the OD pair flow value threshold that does not participate in the AMOEBA calculation

Third, we adopt Monte Carlo simulation based on conditional permutations to evaluate the statistical significance. The default number is 1000.

Finally, after execute the code, export the results to the path you want to save.
outputFile = open('yourpath/file_name.txt','w')
outputFile.write(outputStr)

## Citation 
```
@article{Rantao2023,
author = {Tao, R., Chen, Y., & Thill, J. C.},
doi = {10.1016/j.compenvurbsys.2023.102042},
journal = {Computers, Environment and Urban Systems},
number = {1},
pages = {102042},
title = {{A space-time flow LISA approach for panel flow data}},
url = {https://doi.org/10.1016/j.compenvurbsys.2023.102042},
volume = {106},
year = {2023}}
credits = "Copyright (c) 2018-01 Ran Tao"
license = "New BSD License"
version = "1.0.0"
email = "rtao@usf.edu"
```
