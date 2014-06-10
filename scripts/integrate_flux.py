#!/usr/bin/python

angles = [30, 60]

import sys
import numpy

def main():
    if len(sys.argv) != 2:
        print "Usage:", sys.argv[0], " input_filename"
        exit(1)
    data = []
    with open(sys.argv[1],'r') as f:
        for line in f:
            data.append([float(x) for x in line.split()])
    delta = data[1][0] - data[0][0]
    for i in range(0,len(data) - 1):
        if abs((data[1][0] - data[0][0]) - delta) > 1E-8:
            print "Error: inconsistent delta angles." 
            exit(2)
    totInt = 0
    angInt = 0
    for i in range(0,len(data)):
        totInt += delta * data[i][1]
        angle = data[i][0] * 180 / numpy.pi 
        if angle >= angles[0] and angle <= angles[1]:
            angInt += delta * data[i][1]
            
    #print "Total integral:", totInt, "Between",angles[0],"and",angles[1],angInt 
    print totInt, angInt


if __name__ == "__main__":
    main()
