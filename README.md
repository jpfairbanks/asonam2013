# Case Study using STINGER for social media analysis
In this case study we use Twitter data collected around the impact of hurricane sandy to show how social media can be analyzed using the STINGER graph library.

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
BC can be computed using the commands ...
head BC.csv
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
### From CJ Hutto email
=https://mail.google.com/mail/ca/u/0/#apps/cj/13c1fb675dd7c042=

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

