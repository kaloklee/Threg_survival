# Threg_survival

In this respostory, I use threshold regression for survival Analysis with right-censored data (Lee and Whitmore 2006).

This method attracts me because of its ability to handle non-proportional hazard.  

There is a R package, _threg_, to run this method (Xiao et al. 2015) and I highly recommend it.

In this repository, I coded it up using Stan for me to understand it better. 

This is an example of how it can handle non-proportionality:
![alt text](https://github.com/kaloklee/Threg_survival/blob/main/output.jpg?raw=true)

**Reference**

Lee MLT, Whitmore GA (2006). “Threshold Regression for Survival Analysis: Modeling
Event Times by a Stochastic Process.” Statistical Science, 21(4), 501–513.

Xiao, T., Whitmore, G. A., He, X., & Lee, M.-L. T. (2015). The R Package threg to Implement Threshold Regression Models. Journal of Statistical Software, 66(8), 1–16.





