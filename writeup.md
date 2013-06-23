# Case Study using STINGER for social media analysis

In this case study we use Twitter data collected around the impact of hurricane
sandy to show how social media can be analyzed using the STINGER graph library.
The analysis ranges from the simple degree distribution estimation to tracking
the influential users over time.

We propose a method of analysis for streaming graphs. This method starts from
the graph and leads into plots generated for visualizing the results of the
computational analysis by a human analyst. Every kernel produces either a global
number (global clustering coefficient) or assigns a local number to each vertex
(degree, betweenness centrality). In order to be shown on a plot, this
information must be condensed to two dimensions. One meaningful way to condense
this information is to use a histogram to show the distribution of values. In
the static case this is easy because we can put the sorted distinct kernel
values (bins) on the $X$ axis and the number of vertices with that value the $Y$
axis. This allows an analyst to understand the distribution of the kernel
numbers and thus chose an algorithm to find interesting vertices. This can be
a task such as searching for outliers or segmenting the population. The
knowledge of the distribution is important for the success of many statistical
techniques for these analysis tasks.

In the streaming case we also have a time component that makes this strategy
less viable. We now have a value for each vertex for all time steps. This cannot
be visualized by displaying a single histogram. Thus we propose to take summary
statistics from the distributions and plot them changing over time. This will
allow an analyst to visually inspect the data for significant temporal events.
We can also use time series analysis techniques to find significant events in
these summary statistics.

## Getting Started
+ Install ADAMS-STINGER from Rob.
    - `make install`...
+ Prepare your workspace.
    We assume you have a layout of
    - `code/stinger`
    - `data/`
    - `output`
    - `scripts/`
+ Write config file see [configuration syntax](#config_syntax)
 `code/stinger/config.h`

Here is an example of executing the clustering and connected component kernel
This should build the graph and give you a summary of it to stdout.

TODO: this has changed so you have to cat graph.csv|sandy.c and it will run
everything. This interface will change again when rob makes his function
registration system.

~~~~~~~~~~
>> code/stinger/main data/ output/ -a c
Up and running. Config...
        reading-from: data/
        writing-to: output/
        client-server: stand-alone
        graph-tables-enabled:
                GTYPE_EMAIL_EVENT
                GTYPE_FILE_EVENT
                GTYPE_IM_EVENT
                GTYPE_LOGON
                GTYPE_PRINTER_EVENT
                GTYPE_PROCESS_EVENT
                GTYPE_REGISTRY_EVENT
                GTYPE_URL_EVENT
                GTYPE_USER_EMAIL
                GTYPE_USER_IM
        algorithms-enabled:
                ATYPE_GTRI2
Building Graph...
        Read 1238109 emailEvent rows in 1.705518 seconds
        tree to stinger 1.460194 seconds
Graph(s) created. Running stats...
        Vertices: 662575
        Edges: 2189965
        Done. 0.082313 seconds
Beginning clustering...
        cluster_time 1.017081 seconds
        clusters 2313
        clusters_edges 5528
        clustering_iterations 7
Beginning connected components...
        Done. 0.125424 seconds 49204 components
Algorithms have completed. Closing.
~~~~~~~~~~~~

The -a c is instructing the program to run the **a**lgorithm named
**c**lustering-components. Without the algorithm the program will attempt all
kernels. If the necessary data files are missing, then the program will fail on
initializing them. We can use clustering-components because all of the necessary
data is in the graph specification file that is included by the configuration
file.



## Description of DATA collected
During the event of hurricane Sandy, we collected tweets using the hashtag list
 - sandy
 - huricanesandy
 - frankenstorm
 - zoneA

We collected <#> tweets from <#> distinct users. These were collected
as json objects and then a graph was extracted. The graph of interest
is the at (@) mentions network where two users are connected if one
has mentioned the other in a tweet. We collected a symmetric graph by
adding both directions of the edge. This resulted in a graph
containing <#> vertices and <#> edges. This corresponds to scale <#> and
edgefactor <#>.

## Kernels

### Degree distribution

Describe the commands for computing degree distribution

The static form generates a histogram of the degrees.
For all degrees, the number of vertices with that degree.
The streaming form is an estimate of the log normal parameters, $\mu \sigma$
which assumes that $\log(d_i)$ is normally distributed. This is a simple heavy
tailed distribution and will not be a good fit for a road network.

TODO: insert a degree distribution plot and the logmean log std deviation over
time plots

The vertices are printed if their `VTYPE` is not `VTYPE_NONE`.
This is important if you ever use the `VTYPE_NONE` as a valid type.

### Connected Components

List commands
We found that the distribution of component sizes looks like...

For all component sizes, number of components with that size.

The static form is a histogram of component sizes.
### Betwenness Centrality BC

BC can be computed using the command

   `./code/stinger/main data/ output -a b`


TODO: include `head BC.csv`

We find that the most central users are
the media outlets and government officials. This confirms the findings
in Ediger et al. The interesting things are the dynamics of this list
which will be discussed later <pointer to discussion>.

### community detection
List commands and files generated.
Lets crack open a medium sized community and see which users it has.

## Data Analysis

See the paper published at ASONAM to see what we actually did.

The main goals were to describe the distributions of graph kernel outputs and visualize them.
A successful outlier detection was developed that used both temporal information and the changes in topology.

# Appendix

## config_syntax
