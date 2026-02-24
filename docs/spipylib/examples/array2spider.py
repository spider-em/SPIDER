import numpy
from Spider.Spiderarray import array2spider

size = 65

I = numpy.zeros((size,size))   # create a blank array

for j in range(size):
    for i in range(size):
        I[j][i] = i              # Array[row][column]

array2spider(I, 'tmp001.dat')
