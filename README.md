<!-- README.md is generated from README.Rmd. Please edit that file -->
econet
======

The R package *econet* provides methods for estimating parameter-dependent network centrality measures with linear-in-means models. Both non linear least squares and maximum likelihood estimators are implemented. The methods allow for both link and node heterogeneity in network effects, endogenous network formation and the presence of unconnected nodes. The routines also compare the explanatory power of parameter-dependent network centrality measures with those of standard measures of network centrality. Benefits and features of the econet package are illustrated in the vignette of the package.

Implementing Econet
-------------------

*econet* provides four functions. The first one is *net\_dep*, which allows one to estimate a model of social interactions and compute the relative weighted Katz-Bonacich centralities of the agents. Different behavioral models can be chosen (see section 4 of the vignette for details). Moreover, the hypothesis of homogenous or heterogenous spillovers can be tested. The second function is *boot*, which is built to obtain valid inference when the NLLS estimator with Heckman correction is used. The third function is *horse\_race*, which allows one to compare the explanatory power of parameter-dependent centralities relative to other centrality measures. The forth function is *quantify*, and it is used to assess the effect of control variables in the framework designed by BLP.

The package has at least four merits.

First, it complements the R packages implementing traditional centrality measures for binary networks, *igraph* and *sna*, and weighted networks, *tnet*, by introducing new eigensolutions-based techniques to rank agents' centrality. Second, whereas previous packages, such as *btergm*, *hergm*, the *statnet* suite, and *xergm*, created environments for modeling the statistical processes underlying network formation, *econet* provides the first framework to investigate the socio-economic processes operating on networks (i.e. peer effects). Third, it completes the collection of functions for modeling spatial dependence in cross-sectional data provided by *spdep* and *splm*, by allowing the users to: i) consider the presence of unconnected nodes, and ii) address network endogeneity. Finally, it equips the R archive with routines still unavailable in other commonly used software for the investigation of relational data, such as Matlab, Pajek, Python and Stata.

The examples we use to showcase the functionality of *econet* are contained in the vignette.
