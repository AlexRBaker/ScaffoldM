# ScaffoldM

A coverage-based scaffolder which uses a probabilistic distance measurement () to compare contigs as well as a SSPACE-like (, 2010) link comparison.

## Overview
ScaffoldM emulates features of SSPACE while adding a comparison based on the coverage of contigs. So, first, possible links are evaluated
and all contig pairs with less then k (user defined, deafult=5) links are excluded from further analysis. Then, a random seed contig pair is chosen and 
is extended. This is done by considering each end of the pair and finding other contigs pairs which share an edge. Then, from these plausible pairs an extra decision process
is made. The ratio of the tentative contig pair with the most supporting reads on the number of links for other pairs 
is used to reject making a connection if any of the ratios is >(threshold). After this, the 'sequence space' is searched. The insert size is the size of this space.
The portion of the contig in this space (min(contiglength, insertsize-gapsize)) is then divided by the number of linking reads. These sequence space per link measurements are again
compared by ratios, if any are greater than the threshold no extension is made. Finally, the coverage based comparisons come. The probabilistic distance is calculated as specified by
metabat as the $\sim_\limits_{i=0}^{\infty}|\Psi(i;\mu=\mu_{1},\sigma=\sigma_{1})-\Psi(i;\mu=\mu_{2},\sigma=\sigma_{2})|$. If the distance is >0.95 the connection is rejected as the
coverage is too different between contigs. If the distance is <0.05 then a connect is enforced.

## Installation

Coming soon...

## Example usage

## Help

If you experience any problems using scaffoldm, open an [issue](https://github.com/ecogenomics/scaffoldm/issues) on GitHub and tell us about it.

## Licence and referencing

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/ecogenomics/scaffoldm

This software is currently unpublished

## Copyright

Copyright (c) 2015 Alexander Baker and Michael Imelfort. See LICENSE.txt for further details.
