# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 13:17:14 2020

@author: Qilin Yu
"""
class Point:
    def __init__(self, x, y): 
        self.x = x 
        self.y = y 
        
def ccw(A,B,C):
    return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)

def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)


A = Point(5.461214441214744, 40.81270632864131)
B = Point(5.22585807, 39.7809082)
C = Point(19.999999999999773, 40.0000045)
D = Point(0,40.0000025)

intersect(A,B,C,D)
