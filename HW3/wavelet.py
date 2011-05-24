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

2) Implement Wavelet-Based data fusion method (figure 8 of the Wavelet based fusion methods handout posted on General Resources )

3) Test it on MODIS pre-processed bands that are posted on: http://glasslab.org/u/georgeb/cap2011/
Use MODIS pixel registration chart shown in Figure 1 of MODIS True Color handout (posted on General Resources). 

Note: for LCC, you can use the sample correlation coefficients formula from 
http://en.wikipedia.org/wiki/Correlation_and_dependence
"""


    
    

def _get4widths(width): 
    """
    Professor provided function to get widths
    """
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
                                                                                       
    if s != 0:                                                                                                                 
        south = img[-s:,:]                                                                                                     
        img = N.row_stack((img, south[::-1,:]))                                                                                
        
    if e != 0:                                                                                                                 
        east = img[:,-e:]                                                                                                      
        img = N.column_stack((img, east[:,::-1]))                                                                              
        
    if w != 0:                                                                                                                 
        west = img[:,:w]                                                                                                       
        img = N.column_stack((west[:,::-1], img))                                                                              
        
        
        
        return img
    
    
def enlarge(arr):    
    #ENLARGING ARRAY
    #j = N.array([[1,2],[3,4]])   
    l =N.repeat(arr, 2, axis = 0)
    enlarged =N.repeat(l, 2, axis = 1)
    return enlarged
def Pyramid():    
    #ARRAY T
    x = 4
    y = 4
    len = x* y
    A=N.arange(len).reshape((x,y))
    #print "this is array " ,"\n", A

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
    #print "This is g" , g

                                                     
def UWA_left_split(arr):                                                        
    """
    get left split of array
    """
    col = N.shape(arr)[1]
    arr = arr[:, : col/2]                                                         
    return arr
def UWA_right_split(arr):
    """
    get right split of array
    """
    col = N.shape(arr)[1]
    arr = arr[: , col/2:]    
    return arr
#Beggining of image box
def UWA(arr):
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
            e = N.append(e,(array_slice[1]-array_slice[0])/2.)
        i = i+1
    #This reshapes Jh
    c = N.reshape(c, (rows, columns/2))
    #print "This is c" , c
    #This reshapes Jg
    e = N.reshape(e, (rows,columns/2))
    #print "this is e" , e
    #These two lines input  Jh & Jg into array
    d[:,0:columns/2] = c
    d[:,columns/2: ] = e
    return d
    
 
def UWA2(arr):
    """
    This function takes in the array from UWA and performs Undecimated 
    Wavelet Analysis one more time and returns the left(LL and LH) or right
    (HL and HH) depends on the one needed when function call.
    """
    #arr = N.arange(8).reshape(4,2)
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
            #This takes care of Jhh(LL) or Jgh(HL)
            c = N.append(c,N.average(array_slice))
            #This takes care of Jhg(LH) or Jgg(HH)
            e = N.append(e,array_slice[1]-array_slice[0])
            #print "this is e" , e
        i = i+1
    #This reshapes Jhh(LL) or Jgh(HL)
    c = N.reshape(c, (rows/2, columns))
    #This reshapes Jhg(LH) or Jgg(HH)
    e = N.reshape(e, (rows/2,columns))
    #These two lines input  Jhh & Jhg into array
    d[:rows/2,:] = c
    d[rows/2: ,:] = e
    return d 
    


def window(PAN, MS):
    """
    3X3 sliding window pass in image and a tuple of the window size e.g. (3,3)
    """
    new_arr  = N.array([])
    #im = N.arange(36).reshape(6,6)
    wsize = (3,3)
    dx, dy = wsize
    nx = PAN.shape[1] -dx+1
    ny = PAN.shape[0]- dy+1
    nx = MS.shape[1] -dx+1
    ny = MS.shape[0]- dy+1
    #results= N.array([])
    for i in xrange(ny):
        for j in xrange(nx):
            PANW = PAN[i:i+dy, j:j+dx]
            #PANW = N.ravel(PAN)
            MSW=  MS[i:i+dy, j:j+dx]
            #MSW = N.ravel(MS)
            #print "This is PANW" , PANW , "\n", "This is MSW" , MSW
            LCC = lcc(PANW, MSW)
            
            new_arr = N.append(new_arr,LCC)
    return new_arr

"""Window for LG"""
def windowDF(PAN, MS, HL_PAN, HL_MS):
    """
    3X3 sliding window pass in image and a tuple of the window size e.g. (3,3)
    """
    new_arr  = N.array([])
    #im = N.arange(36).reshape(6,6)
    wsize = (3,3)
    dx, dy = wsize
    nx = PAN.shape[1] -dx+1
    ny = PAN.shape[0]- dy+1
    #nx = MS.shape[1] -dx+1
    #ny = MS.shape[0]- dy+1
    HLA = 0
    #results= N.array([])
    for i in xrange(ny):
        for j in xrange(nx):
            PANW = PAN[i:i+dy, j:j+dx]
            MSW=  MS[i:i+dy, j:j+dx]
            HL_PANW = HL_PAN[i:i+dy, j:j+dx]
            HL_MSW = HL_MS[i:i+dy, j:j+dx]
            print "Thi is HL_PAN" , HL_PANW
            print "This is HL_MS " , HL_MSW
            HL_PANW = N.ravel(HL_PAN)
            HL_MSW = N.ravel(HL_MS)
            LCC = lcc(PANW, MSW)
            print "This is LCC" , LCC
            LG = lg(PANW, MSW)
            print "this is LG" ,LG
            """If LCC is greater than threshold(0.6)"""
            if LCC<0.6:
                LG = 0
            else:
                HLA = LG*HL_PANW[4] +HL_MSW[4] 
            new_arr = N.append(new_arr,HLA)
    return new_arr

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
        #print "This is i PANvar" , PANvar                                             
        #print "This is i MSvar" , MSvar                                               
        LG = MSvar/float(PANvar)                                                      
    return LG   

def lcc(PAN, MS):
    """
    Local Correlation Coefficient function
    """
    #new_arr = N.array([])
    PAN = N.ravel(PAN)
    MS = N.ravel(MS)
    #MS = MS[::-1]
    numerator = 0
    denominator = 0
    for i in range(9):
        numerator += ((PAN[i] - N.average(PAN)) * (MS[i]-N.average(MS)))
        denominator += ((PAN[i]-N.average(PAN))**2 )*((MS[i]-N.average(MS))**2)
    
    #print "numerator" , numerator, "denominator" ,denominator
    result = numerator/float(sqrt(denominator))
    return result


def top_split(arr):
    """
    This is a function to split an array in half along y axis and return top piece
    """
    row = N.shape(arr)[0]
    col= N.shape(arr)[1]
    a= arr[:row/2, :]
    return a
def bottom_split(arr):
    """
    This is a function to split an array in half along y axis and return bottom piece
    """
    row = N.shape(arr)[0]
    col = N.shape(arr)[1]
    b = arr[row/2: , :]
    return b

def main():
    ext = '.raw' if scene == 'coast' else '.dat'
    #PAN = N.arange(144).reshape(12,12)
    PAN = N.fromfile(os.path.join(sourcedir, scene, 'b8'+ext), 'uint8').reshape((M,M))
    MS_R = N.fromfile(os.path.join(sourcedir, scene, 'b3'+ext), 'uint8').reshape((K,K))
    print "This is PAN" , PAN
    print "THis is MS_R" ,MS_R
    sys.exit(0)
    """
    PAN Image functions
    """
    g =UWA(PAN)
    #print "this is g" , g
    #m = N.arange(16).reshape(4,4)
    left_split = UWA_left_split(g)
    #print "this is left split", left_split
    right_split = UWA_right_split(g)
    #print "This is right split " , right_split
    lft_splitted =UWA2(left_split)
    #print "This is lft splitted" , lft_splitted
    rgt_splitted = UWA2(right_split)
    #print "This is right splitted", rgt_splitted
    """
    Lets take lft_splitted and split it into LL and LH
    """
    
    LL_pan = top_split(lft_splitted)
    #print "This is LLpan" , LL_pan
    LH_pan = bottom_split(lft_splitted)
    """
    Lets take rgt_splitted and split it into HL and HH
    """
    HL_pan = top_split(rgt_splitted)
    HH_pan = bottom_split(rgt_splitted)
    
    """
    END of PAN UWA
    """
    
    """
    MS Functions
    """
    """
    MS Red band
    """
    #MS_R = N.arange(36)
    #MS_R = MS_R[::-1].reshape(6,6)
    MS_R = enlarge(MS_R)
    h = UWA(MS_R)
    
    left_split = UWA_left_split(MS_R)
    #print "this is left split", left_split
    right_split = UWA_right_split(MS_R)
    #print "This is right split " , right_split
    lft_splitted =UWA2(left_split)
    #print "This is lft splitted" , lft_splitted
    rgt_splitted = UWA2(right_split)
    #print "This is right splitted", rgt_splitted
    """
    Lets take lft_splitted and split it into LL and LH
    """
    LL_MS_R = top_split(lft_splitted)
    LH_MS_R = bottom_split(lft_splitted)
    """
    Lets take rgt_splitted and split it into HL and HH
    """
    HL_MS_R = top_split(rgt_splitted)
    HH_MS_R = bottom_split(rgt_splitted)
    """
    END of MS Red BAND
    """
    
    """
    MS Green band
    """
    MS_G = N.arange(36)
    MS_G = MS_G[::-1].reshape(6,6)
    MS_G = enlarge(MS_G)
    h = UWA(MS_G)
    
    left_split = UWA_left_split(MS_G)
    #print "this is left split", left_split
    right_split = UWA_right_split(MS_G)
    #print "This is right split " , right_split
    lft_splitted =UWA2(left_split)
    #print "This is lft splitted" , lft_splitted
    rgt_splitted = UWA2(right_split)
    #print "This is right splitted", rgt_splitted
    """
    Lets take lft_splitted and split it into LL and LH
    """
    LL_MS_G = top_split(lft_splitted)
    LH_MS_G = bottom_split(lft_splitted)
    """
    Lets take rgt_splitted and split it into HL and HH
    """
    HL_MS_G = top_split(rgt_splitted)
    HH_MS_G = bottom_split(rgt_splitted)
    """
    END of MS Green BAND
    """
    
    """
    MS Blue band
    """
    MS_B = N.arange(36)
    MS_B = MS_B[::-1].reshape(6,6)
    MS_B = enlarge(MS_B)
    h = UWA(MS_B)
    
    left_split = UWA_left_split(MS_B)
    #print "this is left split", left_split
    right_split = UWA_right_split(MS_R)
    #print "This is right split " , right_split
    lft_splitted =UWA2(left_split)
    #print "This is lft splitted" , lft_splitted
    rgt_splitted = UWA2(right_split)
    #print "This is right splitted", rgt_splitted
    """
    Lets take lft_splitted and split it into LL and LH
    """
    LL_MS_B = top_split(lft_splitted)
    LH_MS_B = bottom_split(lft_splitted)
    """
    Lets take rgt_splitted and split it into HL and HH
    """
    HL_MS_B = top_split(rgt_splitted)
    HH_MS_B = bottom_split(rgt_splitted)
    
    """
    END of MS Blue BAND
    """
    """
    LCC
    """
    """PADDED PAN Image"""
    pad_LL_pan = mirrorpad(LL_pan,(1,1,1,1))
    pad_LH_pan = mirrorpad(LH_pan,(1,1,1,1))
    pad_HL_pan = mirrorpad(HL_pan,(1,1,1,1))
    pad_HH_pan = mirrorpad(HH_pan,(1,1,1,1))
    
    """PADDED MS Red Band"""
    pad_LL_MS_R = mirrorpad(LL_MS_R, (1,1,1,1)) 
    pad_LH_MS_R = mirrorpad(LH_MS_R, (1,1,1,1)) 
    pad_HL_MS_R = mirrorpad(HL_MS_R, (1,1,1,1)) 
    pad_HH_MS_R = mirrorpad(HH_MS_R, (1,1,1,1)) 
    
    """PADDED MS Green Band"""
    pad_LL_MS_G = mirrorpad(LL_MS_G, (1,1,1,1)) 
    pad_LH_MS_G = mirrorpad(LH_MS_G, (1,1,1,1)) 
    pad_HL_MS_G = mirrorpad(HL_MS_G, (1,1,1,1)) 
    pad_HH_MS_G = mirrorpad(HH_MS_G, (1,1,1,1)) 
    
    """PADDED MS Blue Band"""
    pad_LL_MS_B = mirrorpad(LL_MS_B, (1,1,1,1)) 
    pad_LH_MS_B = mirrorpad(LH_MS_B, (1,1,1,1)) 
    pad_HL_MS_B = mirrorpad(HL_MS_B, (1,1,1,1)) 
    pad_HH_MS_B = mirrorpad(HH_MS_B, (1,1,1,1)) 
    
    """LCC LL PAN with LL MS Red band"""
    LCC_LL_MS_R = window(pad_LL_pan,pad_LL_MS_R)
    """LCC LL PAN with LL MS Green band"""
    LCC_LL_MS_G = window(pad_LL_pan,pad_LL_MS_G)
    """LCC LL PAN with LL MS Blue band"""
    LCC_LL_MS_B = window(pad_LL_pan,pad_LL_MS_B)
    
    
    
    HLRED = windowDF(pad_LL_pan, pad_LL_MS_R,pad_HL_pan,pad_HL_MS_R)
    HHRED = windowDF(pad_LL_pan,pad_LL_MS_R,pad_HH_pan,pad_HH_MS_R)
    LHRED = windowDF(pad_LL_pan,pad_LL_MS_R,pad_LH_pan,pad_LH_MS_R)
    
    HLGREEN = windowDF(pad_LL_pan, pad_LL_MS_G, pad_HL_pan,pad_HL_MS_R)
    HHGREEN = windowDF(pad_LL_pan,pad_LL_MS_G,pad_HH_pan,pad_HH_MS_G)
    LHGREEN = windowDF(pad_LL_pan,pad_LL_MS_G,pad_LH_pan,pad_LH_MS_G)
    
    
    HLBLUE = windowDF(pad_LL_pan, pad_LL_MS_B, pad_HL_pan,pad_HL_MS_B)
    HHBLUE = windowDF(pad_LL_pan,pad_LL_MS_B,pad_HH_pan,pad_HH_MS_B)
    LHBLUE = windowDF(pad_LL_pan,pad_LL_MS_B,pad_LH_pan,pad_LH_MS_B)
    
    
    
    
    
    """
    END of LCC
    """
    
    """
    
    """
if __name__== "__main__":
    main()
    #A=N.arange(16).reshape((4,4))
    #print (neighbors(A,0,0))
    #Pyramid()
    #PAN = N.arange(144).reshape(12,12)
    #MS = N.arange(144)
    #MS = MS[::-1].reshape(12,12)
    #UWA(PAN)
    #UWA2()
    #r = N.arange(64).reshape(8,8)
    #u = N.arange(64).reshape(8,8)
    
    #window(r, u)
    #lcc()