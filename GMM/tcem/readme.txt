Matlab code for the algorithm of the following paper

G.Lee and C.Scott, "EM algorithms for multivariate Gaussian mixture models with truncated and censored data"

Author : Gyemin Lee
Date : August, 2010.

This software is inteded for noncommerial purposes.

----------------------------------------

--Demo
run 'demo.m'



--Files

demo.m : demo on a bivariate three-components Gaussian mixture sample


em_censor_pattern.m : find the censored data pattern


em_std.m	: standard EM
em_std_post.m


em_tc1.m	: truncated and censored data EM 
		  currently implementation considered the truncation is only on the first coordinate
em_tc1_post.m
em_censor_post.m


tmvn_m3.m : mean and covariance of a truncated multivariate Gaussian distribution


init_kmeans.m : EM algorithm is initialized from k-means


util_*.m : helper functions
