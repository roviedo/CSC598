#!/usr/bin/python

import numpy as N
from math import sqrt 
"""
1) Implement Pyramid-Based fusion method (figure 9 of the Wavelet based fusion methods handout posted on General Resources )

2) Implement Wavelet-Based data fusion method (figure 8 of the Wavelet based fusion methods handout posted on General Resources )

3) Test it on MODIS pre-processed bands that are posted on: http://glasslab.org/u/georgeb/cap2011/
Use MODIS pixel registration chart shown in Figure 1 of MODIS True Color handout (posted on General Resources). 

Note: for LCC, you can use the sample correlation coefficients formula from 
http://en.wikipedia.org/wiki/Correlation_and_dependence
"""

def Pyramid():
    
    #RAVEL EXAMPLE
    X = N.matrix('1 2 3 4; 5 6 7 8;9 10 11 12;13 14 15 16')
    flat = N.ravel(X)
    print "This is the flatten array" , flat
    
    #RESHAPE EXAMPLE
    B=N.arange(48).reshape(4,4,3)
    flat2 = N.ravel(B)
    #print "This is B array ", "\n" , B 
    #print "This is flatten 4 x 4 x 3:" , "\n" , flat2 
    #C= N.mat(B.copy())
    #print "THis is matrix" , C
    """
    """
    C = N.arange(4).reshape(2,2)
    a = N.array([[1, 2],[3, 4]])  
    
    #REPEAT EXAMPLE
    #b = N.zeros((4,4))  
    #b[:,1:-1] = N.repeat(a, 2,axis=0)
    #print "This is b" , b
    #c = N.arange(4)
    """
    c = N.array([2 ,3, 5, 6])
    d = N.array([])
    for i in c:
        #d = c.repeat(i)
         n = N.repeat(i , 2)
         d= N.append(d, n)
    print "This is c \n" , d
    d = d.reshape(2, 4)
    print "this is d ", "\n", d
    """
    e = N.array([[1,2],[3,4]])
    f = N.zeros((4, 4))
    #f[:,1 : -1] = N.repeat(e, 2, axis=0)

    #ENLARGING ARRAY
    j = N.array([[1,2],[3,4]])   
    l =N.repeat(j, 2, axis = 0)
    k =N.repeat(l, 2, axis = 1)
    
    #ARRAY T
    x = 4
    y = 4
    len = x* y
    A=N.arange(len).reshape((x,y))
    print "this is array " ,"\n", A

    #PAD ARRAY WITH ZEROS
    C = N.insert (A, (0,x), 0, axis=1)
    C = N.insert (C, (0,y), 0, axis=0)
    #print "This is array C " ,"\n", C
    
    
    
    
    #REDUCING ARRAY 
    g = N.array([])
    for i in range(x):
        for j in range(y):
            if i%2 == 0 and j%2==0:
                slice_array =A[i:(2+i) , j:(2+j)]
                g = N.append(g, N.average(slice_array))
    print "This is g" , g


def UWA(arr):
    """
    Beggining of image box
    """
    #arr = N.arange(16).reshape(4,4)
    rows = N.shape(arr)[0]
    columns = N.shape(arr)[1]
    c = N.array([])
    e = N.array([])
    d = N.zeros([rows,columns])
    B = N.ravel(arr)
    i = 0
    while (i<len(B)):
        if i%2 == 0:
            array_slice = B[i:(2+i)]
            #This takes care of Jh
            c = N.append(c,N.average(array_slice))
            #This takes care of Jg
            e = N.append(e,(array_slice[1]-array_slice[0])/2)
        i = i+1
    #This reshapes Jh
    c = N.reshape(c, (rows, columns/2))
    #This reshapes Jg
    e = N.reshape(e, (rows,columns/2))
    #These two lines input  Jh & Jg into array
    d[:,0:columns/2] = c
    d[:,2: ] = e
    return d

def UWA_left_split(arr):
    """
    get left split of array
    """
    arr = a[:, : col/2]
def UWA_right_split(arr):
    arr = a[: , col/2:]
    
    
def UWA2(arr):
    """
    pass in left half and right half seperately
    """
    arr = N.arange(8).reshape(4,2)
    rows = N.shape(arr)[0]
    columns = N.shape(arr)[1]
    c = N.array([])
    e = N.array([])
    d = N.zeros([rows,columns])
    B = N.ravel(arr)
    i = 0
    while (i<len(B)):
        if i%2 == 0:
            array_slice = B[i:(2+i)]
            #This takes care of Jhh
            c = N.append(c,N.average(array_slice))
            #This takes care of Jhg
            e = N.append(e,array_slice[1]-array_slice[0])
        i = i+1
    #This reshapes Jhh
    c = N.reshape(c, (rows/2, columns))
    #This reshapes Jhg
    e = N.reshape(e, (rows/2,columns))
    #These two lines input  Jhh & Jhg into array
    d[:rows/2,:] = c
    d[rows/2: ,:] = e
    return d 
    

#3X3 sliding window pass in image and a tuple of the window size e.g. (3,3)
def window():
    im = N.arange(36).reshape(6,6)
    wsize = (3,3)
    dx, dy = wsize
    nx = im.shape[1] -dx+1
    ny = im.shape[0]- dy+1
    results= N.array([])
    for i in xrange(ny):
        for j in xrange(nx):
            g =  im[i:i+dy, j:j+dx]
            print "this is window" , g
            g = N.ravel(g)
    print "this is g " , g
def lcc():
    new_arr = N.array([])
    PAN = N.arange(9)
    MS = N.arange(9)
    MS = MS[::-1]
    numerator = 0
    denominator = 0
    for i in PAN:
        numerator += ((PAN[i] - N.average(PAN)) * (MS[i]-N.average(MS)))
        denominator += ((PAN[i]-N.average(PAN))**2 )*((MS[i]-N.average(MS))**2)
    
    print "numerator" , numerator, "denominator" ,denominator
    result = numerator/float(sqrt(denominator))
    print "this is result" , result
    
if __name__== "__main__":
    #A=N.arange(16).reshape((4,4))
    #print (neighbors(A,0,0))
    #Pyramid()
    #UWA(A)
    #UWA2()
    window()
    #lcc()