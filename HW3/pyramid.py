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

def enlarge(arr):
    """
    ENLARGING ARRAY to twice size                                                                
    """
    #j = N.array([[1,2],[3,4]])                                                     
    l =N.repeat(arr, 2, axis = 0)                                                   
    enlarged =N.repeat(l, 2, axis = 1)                                              
    return enlarged  

def reduce(A):
    """
    Reduce image to half size
    """
    x = N.shape(A)[0]
    y = N.shape(A)[1]
    g = N.array([])
    for i in range(x):
        for j in range(y):
            if i%2 == 0 and j%2==0:
                slice_array =A[i:(2+i) , j:(2+j)]
                g = N.append(g, N.average(slice_array))
                
    gg = g.reshape(3,3)            
    return gg

def _get4widths(width):                                                         
    try:                                                                        
        width = int(width)                                                      
        width = (width, width, width, width)                                    
    except TypeError:                                                           
        width = tuple(width)                                                    
        if len(width) != 4:                                                     
            raise ValueError("width must be either an integer or a 4-tuple")    
        if any([x < 0 for x in width]):                                         
            raise ValueError("negative value in width=%s" % width)              
        return width 
    
def mirrorpad(img, width):                                                                                                       
    """Return the image resulting from padding width amount of pixels on each                                                    
    sides of the image img.  The padded values are mirror image with respect to                                                  
    the borders of img.  Width can be an integer or a tuple of the form (north,                                                  
    south, east, west).                                                                                                          
    """                                                                                                                          
    n, s, e, w = _get4widths(width)                                                                                              
    row = N.shape(img)[0]                                                                                                        
    
    col = N.shape(img)[1]                                                                                                        
    
    if n != 0:                                                                                                                   
        north = img[:n,:]                                                                                                        
        img = N.row_stack((north[::-1,:], img))                                                                                  
        print "This is img north", img                                                                                           
    if s != 0:                                                                                                                   
        south = img[-s:,:]                                                                                                       
        img = N.row_stack((img, south[::-1,:]))                                                                                  
        print "This is img south", img                                                                                           
    if e != 0:                                                                                                                   
        east = img[:,-e:]                                                                                                        
        img = N.column_stack((img, east[:,::-1]))                                                                                
        print "This is img east", img                                                                                            
    if w != 0:                                                                                                                   
        west = img[:,:w]                                                                                                         
        img = N.column_stack((west[:,::-1], img))                                                                                
        
    return img          
    
#3X3 sliding window pass in image and a tuple of the window size e.g. (3,3)
def window():
    im  = N.arange(36).reshape(6,6)
    im2  = N.arange(36).reshape(6,6)
    wsize = (3,3)
    dx, dy = wsize
    nx = im.shape[1] -dx+1
    ny = im.shape[0]- dy+1
    
    results= N.array([])
    for i in xrange(ny):
        for j in xrange(nx):
            g =  im[i:i+dy, j:j+dx]
            h =  im2[i:i+dy, j:j+dx]
            #PAN = N.array([11,23,45,65,78,22,89,76,90])
            #MS = N.array([56,54,78,34,69,21,48,98,46])
            LG = lg(g,h)
    print "this is LG" , LG

def lcc():
    #new_arr = N.array([])
    PAN = N.arange(9)
    print "PAN" , PAN
    MS = N.arange(9)
    print "MS" , MS
    #MS = MS[::-1]
    numerator = 0
    denominator = 0
    for i in PAN:
        numerator += ((PAN[i] - N.average(PAN)) * (MS[i]-N.average(MS)))
        denominator += ((PAN[i]-N.average(PAN))**2 )*((MS[i]-N.average(MS))**2)
    
    print "numerator" , numerator, "denominator" ,denominator
    result = numerator/float(sqrt(denominator))
    print "this is result" , result

    
def lg(PAN, MS):
    
    PAN = N.ravel(PAN)
    MS = N.ravel(MS)
    PANavg = N.average(PAN)
    MSavg = N.average(MS)
    PANvar = 0
    MSvar = 0
    for i in range(9):
        PANvar += (PAN[i] - PANavg)
        MSvar += (MS[i] -   MSavg)
        print "This is i PANvar" , PANvar
        print "This is i MSvar" , MSvar
        LG = MSvar/float(PANvar)
    return LG    



def main():
    PAN = N.arange(36).reshape(6,6)
    MS_RED = N.arange(16).reshape(4,4)
    MS_GREEN = N.arange(16).reshape(4,4)
    MS_BLUE = N.arange(16).reshape(4,4)
    red_PAN = reduce(PAN)
    print "red_pan" , red_PAN
    enlarged_PAN =enlarge(red_PAN)
    print "this is enlarged", enlarged
     
    sub_PAN = N.subtract(PAN ,enlarged)
    print "This is sub" ,sub
    
    
    
if __name__== "__main__":
    main()    
    #Pyramid()
    #window()
    ##lcc()
    #lg()