# Case Study using STINGER for social media analysis
In this case study we use Twitter data collected around the impact of hurricane sandy to show how social media can be analyzed using the STINGER graph library.

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
        EBlocks[ETYPE_EMAIL_ALIAS]: 0
        EBlocks[ETYPE_EMAIL_FROM]: 725575
        EBlocks[ETYPE_EMAIL_TO]: 0
        EBlocks[ETYPE_ACCESSED]: 0
        EBlocks[ETYPE_LOGON]: 0
        EBlocks[ETYPE_LOGOFF]: 0
        EBlocks[ETYPE_PRINT]: 0
        EBlocks[ETYPE_IM_FROM]: 0
        EBlocks[ETYPE_IM_TO]: 0
        EBlocks[ETYPE_IM_ALIAS]: 0
        EBlocks[ETYPE_CLIPBOARD]: 0
        EBlocks[ETYPE_REGISTRY]: 0
        Consistency 0
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

The -a c is instructing the program to run the **a**lgorithm named **c**lustering-components. Without the algorithm the program will attempt all kernels. If the necessary data files are missing, then the program will fail on initializing them. We can use clustering-components because all of the necessary data is in the graph specification file that is included by the configuration file.



## Description of DATA collected
During the event of hurricane Sandy, we collected tweets using the hashtag list
 - sandy
 - huricanesandy
 - frankenstorm
 - zoneA

We collected <#> tweets from <#> distinct users.
These were collected as json objects and then a graph was extracted.
The graph of interest is the at (@) mentions network where two users are connected if one has mentioned the other in a tweet. We collected a symmetric graph by adding both directions of the edge. This resulted in a graph containing <#> vertices and <#> edges.

## Kernels

### Degree distribution
Describe the commands for computing degree distribution

### Connected Components
List commands 
We found that the distribution of component sizes looks like...

### Betwenness Centrality BC
BC can be computed using the command
    
   `./code/stinger/main data/ output -a b`


TODO: include `head BC.csv`
We find that the most central users are the media outlets and government officials. This confirms the findings in Ediger et al. The interesting things are the dynamics of this list which will be discussed later <pointer to discussion>.

### community detection
List commands and files generated.
Lets crack open a medium sized community and see which users it has.

## Data Analysis

### BC
Measuring the change in permutations can be done with 
Generalized Kendall's Tau
This seems to be the most useful
$\sum_{i<j} p_i p_j [\sigma(i) > \sigma(j)]$
Total number of inversions with positional weights.
James thinks that the positional weights should be exponentially decaying from the top. That way we can ignore the differences at the bottom which will be mostly noise.

### Communities

## Colaboration with CJ Hutto

### Predicting truth using differential centralities
We can compute BC on the large network and the restrict to a rumor (using a regex). Look for a correlation between the difference of 
BC_all, BC_topic_
and the truth or falsity of the rumor.
if the source of deception is outside of twitter allows the first person to be truthfull. Looking at similar tweets at the same time, then there are multiple sources.

### Tracking Spread of topics
Put all of the rumor tweets into a timeline and track the breadth of the timeline and the Min, Q1, Mean Q2, Q3, Max BC of the users tweeting it is. 

Tracking breadth of the filtered network using only edges that go forward in time.

* Filter network using either hashtags or regex
* Do a single edge pass to find the first edges
    + There might be multiple sources of the rumor if it is generate offline
* BFS from there using only edges that increase in time to estimate diameter at each batch step
    + can we measure if it is growing out one vertex at a time?
    + If we don't mark vertices as we see them we can measure if people focus on something for a while or tweet once and move on.
* Keep track of the edges per minute as we add edges.
    + this could be a mean of $t_k-t_{0}$ for each batch of k edges
    + or a sliding window $t_i - t_{i-k}$ for all $e_i$ in the series

### Tracking BC over the course of a rumor spread
We could show that a tweet goes from a medium BC node to a high BC node and then fans out if it is false and spreads. Perhaps if it is true and spreads, then the BC stays low for a while at the begining and then becomes steadily higher until it fans out. This corresponds to falsehoods spreading by decieving an influential source, and then freeloading on his influence. 

A truth will spread by convincing people. Thus we should look at the edge betweeness, if a transfer of information leads to a big disemination and that is the only source.
Then it might be more likely false.

### Examples of Rumors

We can get some ground truth rumors from FEMA.
March madness using UIUC
Lags and Leads of common events.

### From CJ Hutto email

[gmail link](https://mail.google.com/mail/ca/u/0/#apps/cj/13c1fb675dd7c042=)

Something that I think would fit into a collaborative study would be something like comparing various computed measures of influence, such as:

   - number of followers
   - follower-to-following ratio (FFR)
   - lists-to-followers ratio (LFR) http://t.co/TZNgEP3e 
   - number of RTs
   - avg breadth/distance of RTs, 
   - number of @mentions (other than RTs)
   - network centrality (raw)
   - centrality w/in @mention network for given hashtag topic(s)
   - Klout influence score
   - PeerIndex influence score
   - tie strength (as described by Gilbert 2012).

And compare it with a "ground truth" measure of influence... where we ask folks to self report (via survey questionnaire) who in the network was influential both generally and for the given topic. 

I can handle getting IRB approval and design the questionnaire. We can use your set of users in the #Sandy @mention network as our sample (we already have betweeness centrality, so would just need to compute the other influence metrics for comparison). 

We might aim for a publication at IEEE SocialCom conference...

## Conclusions

# Appendix

## config_syntax
