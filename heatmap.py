from __future__ import division

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import re

# Open out.txt for reading
with open('out.txt', 'r') as f:

    # 2d list to store xvals, yvals, and counts
    xvals = [[0 for x in range(32)] for y in range(32)]
    yvals = [[0 for x in range(32)] for y in range(32)]
    count = [[0 for x in range(32)] for y in range(32)]
    
    line = f.readline();
    while line:
    
        xy = map(int, re.findall(r'\d+', line))
        
        if (len(xy) < 4):
            line = f.readline();
            continue
        
        x = xy[0]
        y = xy[1]
        xval = xy[2]
        yval = xy[3]
        
        xvals[x][y] += xval
        yvals[x][y] += yval
        count[x][y] += 1
        
        line = f.readline();
    
    # Calculate average x and y values
    for i in range(32):
        for j in range(32):
            if (count[i][j] != 0):
                xvals[i][j] = xvals[i][j]/count[i][j]
                yvals[i][j] = yvals[i][j]/count[i][j]
    
    # 2d list to store average distances
    distances = [[0 for x in range(32)] for y in range(32)]
    
    # Calculate average distances
    for i in range(32):
        for j in range(32):
            if (count[i][j] != 0):
                x = i
                y = j
                xval = xvals[i][j]
                yval = yvals[i][j]
                distances[i][j] = math.sqrt((xval-x)**2 + (yval-y)**2)
           
    # Find max distance
    max = 0
    for i in range(32):
        for j in range(32): 
            if (count[i][j] != 0):
                if (distances[i][j] > max):
                    max = distances[i][j]
    
    for i in range(32):
        for j in range(32): 
            if (count[i][j] != 0 or distances[i][j] == 0):
                if (distances[i][j] > 4 or distances[i][j] == 0):
                    distances[i][j] = 4
            else:
                distances[i][j] = 1
    
    arr = np.array(distances)
    plt.imshow(arr, origin='lower', cmap='cividis')
    plt.colorbar()
    plt.savefig('plots/heatmap.png')