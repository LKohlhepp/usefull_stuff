# Fancy Counting

This modules makes historgamms, for values with uncertainties. Resulting in a histogramm, in which the counts for each bin also have uncertainties. 

## How to use

Call the function `fancy_bin`.

Its arguments are
- `x`: Are the values that are binned into a histogramm
- `err`: Are the uncertainties of `x` (1 &sigma;)
- `bins`: Are the bin boarders. n+1 boarders need to be supplied for a total of n bins. 

and some koarguments:
- `errfunc`: Is the probablity distribution used to discribe the error. It needs to be the cumulative distribution function of the error, with x to infinity => p = 1, the mean is assumed to be at x = 0 and for sigma = 1. Default: `cdf_gauss`: Gaussianerrorfunction
- `impossiblity_thresehold`: Chances smaller then this for a data point to be in a specific bin are neglegted. This speeds up the calculation, by neglecting very small chances for every point to be in every bin. Default is 10^-4

it returns a numpy array with `size` is `n_bins, 2`. Where in the 0 of the 1st axis is the mean and in the 1 is the error (gaussian) of the counts for the bin. And the 0th axis iterates over all bins. 

## How it works

For each given value the pobabitliy to belong to any bin is calculated. Probabilities below a certain thereshold are set to 0 and are neglected at this point. 
Then for each bin the probablity for any number of archivable counts is determined. Or in other words from 0 counts in this bin to all values that have a non-zero probablity are in this bin, the probablity is determind. 
With these probablities the mean and variance can be calculated.

mean = \<x\> = &Sigma; p(x)x

with x being the number of counts in this bin and p(x) being the probablity for it. 

var = \<x^2\> - \<x\>^2

and 

\<x^2\> = &Sigma; p(x)x^2.
