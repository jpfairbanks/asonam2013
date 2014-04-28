asonam2013
==========

Code accompanying my 2013 ACM ASONAM paper

Getting Started
===============

Each python file in the code directory is a standalone program that takes command line arguments
calling 

      clustering_coefficients.py --help
    
will explain the arguments that it takes including a filename that 
contains data. The code in the stinger directory is the high performance graph code.
Documentation for STINGER can be found at [[stingergraph.com]]. It should install with 

      cd code/stinger
      make all

You should run the stinger code in order to generate the vertex features for the follow on analysis
using the python code. The STINGER code of interest is in the file ~sandy.c~ and can be run with 

      code/stinger/sandy --help
      code/stinger/sandy -n 10 -b 1000 -i 1000 -p ./ < data.csv

Which will describe the options. You must specify the locations of data files and the computations that you 
are interesting in computing, some options are betweenness centrality (approximate) and clustering coefficient
(exact, streaming).

Organization
============

The python code is organized into module files and command line programs.
The following files are the modules than can be imported for future programs.

    - kernelio.py : handles the Input/Output for kernel inputs
    - kernel_analysis.py : holds the functions common to many kernels (vertex features)
    - plotting.py : generic functions for making plots
    - paper_figures.py : specific functions for making plot for the paper presented at asonam.
    
The following programs are intended to be used on the command line. Then can be modified in order to make new programs.
  
    - multivariate.py : code for understanding multiple features at once.
    - anomalies.py : betweenness centrality based outlier detection
    - cc_anomalies.py : clustering coefficient based outlier detection.
    - clustering_coefficients.py : present an analysis of the clustering coefficients.
    - betweenness_centrality.py : present an analysis of the betweenness centrality.
    
