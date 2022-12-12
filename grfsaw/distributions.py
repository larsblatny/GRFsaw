# Author: Lars Blatny
import numpy as np
from scipy.special import gamma
from scipy.stats import rv_continuous, truncnorm

class gamma_distribution(rv_continuous):
    '''
    Gamma probability distribution function used for sampling the magnitude of the wavevector.
    '''
    def _pdf(self, x, b, mean):
        return (b+1)**(b+1) * (x/mean)**b * np.exp( -(b+1)*(x/mean) ) / ( gamma(b+1) * mean )

def normal_distribution(mean, std, lower=0, upper=1e10):
    '''
    Gaussian probability distribution function used for sampling the magnitude of the wavevector. 
    '''
    return truncnorm( (lower-mean)/std, (upper-mean)/std, loc=mean, scale=std )
