#!/usr/bin/python                                                               
import numpy as N                                                               
import os.path                                                                  
import scipy.signal                                                             
import scipy.interpolate                                                        
import matplotlib.pyplot as plt                                                 
import matplotlib.cm as cm                                                      
from math import sqrt                                                           
import sys                                                                      
sourcedir = '/home/raul/CSC598/HW3/'                                            
scene = 'Agriculture'   # 'mountain', 'coast'                                   
K = 1024                                                                        
M = 2048                                                                        
#bands = ['b1', 'b2', 'b3', 'b4', 'b5', 'b7']                                   
bands = ['b3','b2','b1']  


"""
1) Implement Pyramid-Based fusion method (figure 9 of the Wavelet based fusion methods handout posted on General Resources )
   pyramid.py is the file for this method

2) Implement Wavelet-Based data fusion method (figure 8 of the Wavelet based fusion methods handout posted on General Resources )
   wavelet.py is the file for this method

3) Test it on MODIS pre-processed bands that are posted on: http://glasslab.org/u/georgeb/cap2011/
Use MODIS pixel registration chart shown in Figure 1 of MODIS True Color handout (posted on General Resources). 
   Not yet implemented

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
    #print "this is length of A" , x
    y = N.shape(A)[1]
    g = N.array([])
    for i in range(x):
        for j in range(y):
            if i%2 == 0 and j%2==0:
                slice_array =A[i:(2+i) , j:(2+j)]
                g = N.append(g, N.average(slice_array))
                #print "length of g " , len(g)
                
    gg = N.reshape(g,(x/2,y/2))            
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
        #print "This is img north", img                                                                                           
    if s != 0:                                                                                                                   
        south = img[-s:,:]                                                                                                       
        img = N.row_stack((img, south[::-1,:]))                                                                                  
        #print "This is img south", img                                                                                           
    if e != 0:                                                                                                                   
        east = img[:,-e:]                                                                                                        
        img = N.column_stack((img, east[:,::-1]))                                                                                
        #print "This is img east", img                                                                                            
    if w != 0:                                                                                                                   
        west = img[:,:w]                                                                                                         
        img = N.column_stack((west[:,::-1], img))                                                                                
        
    return img          
    

def mirrorunpad(img, width):                                                   
    """
    Return unpadded image of img padded with :func:mirrorpad.               
    """                                                                        
    n, s, e, w = _get4widths(width)                                            
    # index of -0 refers to the first element                                  
    if s == 0:                                                                 
        s = img.shape[0]                                                       
    else:                                                                      
        s = -s                                                                 
    if e == 0:                                                                 
        e = img.shape[1]                                                       
    else:                                                                      
        e = -e                                                                 
    return img[n:s, w:e]  

def window(PAN, MS, sub_PAN):
    """
    3X3 sliding window pass in image and a tuple of the window size e.g. (3,3)
    """
    #PAN  = N.arange(36).reshape(6,6)
    #MS  = N.arange(36).reshape(6,6)
    #sub_PAN = N.arange(36).reshape(6,6)
    length = N.shape(PAN)[0]
    width = N.shape(PAN)[1]
    print "this is length", length, "this is width" ,width
    wsize = (3,3)
    dx, dy = wsize
    nx = PAN.shape[1] -dx+1
    ny = PAN.shape[0]- dy+1
    f = 0
    results= N.array([])
    for i in xrange(ny):
        for j in xrange(nx):
            g =  PAN[i:i+dy, j:j+dx]
            h =  MS[i:i+dy, j:j+dx]
            sub = sub_PAN[i:i+dy,j:j+dx]
            #PAN = N.array([11,23,45,65,78,22,89,76,90])
            #MS = N.array([56,54,78,34,69,21,48,98,46])
            gg = N.ravel(g)
            hh = N.ravel(h)
            subb = N.ravel(sub)
            LCC = lcc(gg,hh)
            LG = lg(gg,hh)
            if LCC<=0.6:
                LG = 0
            else:
                f = LG*subb[4] + hh[4] 
            results = N.append(results, f)
    new_arr = N.reshape(results,(length-2, width-2))
    return new_arr
    

def lcc(PAN, MS):
    """
    Local Correlation Coefficient Function
    """
    #new_arr = N.array([])
    #PAN = N.arange(9)
    #MS = N.arange(9)
    #MS = MS[::-1]
    PAN = N.ravel(PAN)
    MS = N.ravel(MS)
    numerator = 0
    denominator = 0
    for i in range(9):
        numerator += ((PAN[i] - N.average(PAN)) * (MS[i]-N.average(MS)))
        denominator += ((PAN[i]-N.average(PAN))**2 )*((MS[i]-N.average(MS))**2)
    if (float(sqrt(denominator))) == 0:
        result = 0
    else:
        #print "numerator" , numerator, "denominator" ,denominator
        result = numerator/float(sqrt(denominator))
        #print "this is result" , result
    return result
    
def lg(PAN, MS):
    """
    Local gain function
    """
    PAN = N.ravel(PAN)
    MS = N.ravel(MS)
    PANavg = N.average(PAN)
    MSavg = N.average(MS)
    PANvar = 0
    MSvar = 0
    for i in range(9):
        PANvar += (PAN[i] - PANavg)
        MSvar += (MS[i] -   MSavg)
        if (float(PANvar))==0:
            LG = 0
        else:
            LG = MSvar/float(PANvar)
    return LG    



def main():
    """
    Load Landsat files
    """
    ext = '.raw' if scene == 'coast' else '.dat'
    
    #PAN = N.fromfile(os.path.join(sourcedir, scene, 'b8'+ext), 'uint8').reshape((M,M)) 
    #MS_RED = N.fromfile(os.path.join(sourcedir, scene, 'b3'+ext), 'uint8').reshape((K,K))
    #MS_GREEN = N.fromfile(os.path.join(sourcedir, scene, 'b2'+ext), 'uint8').reshape((K,K)) 
    #MS_BLUE = N.fromfile(os.path.join(sourcedir, scene, 'b1'+ext), 'uint8').reshape((K,K))
    
    """
    Testing arrays below
    """
    PAN = N.arange(144).reshape(12,12)
    MS_RED = N.arange(36).reshape(6,6)
    MS_GREEN = N.arange(36).reshape(6,6)
    MS_BLUE = N.arange(36).reshape(6,6)
    """
    Here is our reduced PAN image
    """
    red_PAN = reduce(PAN)
    
    """
    Below we get our enlarged PAN after we reduced it 
    """
    enlarged_PAN =enlarge(red_PAN)
    
    """
    This here is the function that subtract PAN - enlarged PAN to get edges
    """
    sub_PAN = N.subtract(PAN ,enlarged_PAN)
    
    
    """
    These are the three bands enlarged to match size of PAN image
    """
    enlarged_MS_RED = enlarge(MS_RED)
    enlarged_MS_GREEN = enlarge(MS_GREEN)
    enlarged_MS_BLUE = enlarge(MS_BLUE)
    
    """
    These are the padded PAN and MS images
    """
    pad_enlarged_PAN = mirrorpad(enlarged_PAN, (1,1,1,1))
    pad_enlarged_MS_RED = mirrorpad(enlarged_MS_RED, (1,1,1,1))
    pad_enlarged_MS_GREEN = mirrorpad(enlarged_MS_GREEN, (1,1,1,1))
    pad_enlarged_MS_BLUE = mirrorpad(enlarged_MS_BLUE, (1,1,1,1))
    pad_sub_PAN = mirrorpad(sub_PAN , (1,1,1,1))
    
    """
    Functions below are after we multiplied LG with edges and added it 
    to MS image
    """
    
    HR_MS_RED = window( pad_enlarged_PAN, pad_enlarged_MS_RED,pad_sub_PAN)
    print "HiRes_MS_RED Band" , HR_MS_RED
    HR_MS_GREEN = window(pad_enlarged_PAN,pad_enlarged_MS_GREEN,pad_sub_PAN)
    print "HiRes_MS_GREEN Band" , HR_MS_GREEN
    HR_MS_BLUE = window(pad_enlarged_PAN,pad_enlarged_MS_BLUE,pad_sub_PAN)
    print "HiRes_MS_BLUE Band" , HR_MS_BLUE
    
    
    
if __name__== "__main__":
    main()    
    