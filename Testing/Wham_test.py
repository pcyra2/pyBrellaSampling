import numpy

import pyBrellaSampling.wham as wham
import numpy as numpy
def sin(x):
    y = numpy.sinh(x)
    return y

def test_autocorrelation():
    x = [i for i in range(1,100)]
    y = sin(x)
    correlation = wham.autocorrelate(y)
    assert correlation == 1, "Autocorrelation function is not working"
