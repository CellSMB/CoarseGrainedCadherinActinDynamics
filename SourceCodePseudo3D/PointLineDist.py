# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 12:36:49 2020

@author: Qilin Yu
"""
import math

class Point:
    def __init__(self, x, y): 
        self.x = x 
        self.y = y 

def DistSlope(A,B,P):
    slope = (B.y-A.y)/(B.x-A.x)
    intersect = (B.y-B.x*slope)
    dist = abs(-slope*P.x+P.y-intersect)/math.sqrt(slope**2+1)
    return slope, dist

A = Point(1,1)
B = Point(2,2)
P = Point(1,2)

dist = DistSlope(A,B,P)







