import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

# Parse command line args.
xarg = -1
yarg = -1
try:
    xarg = sys.argv[1]
    yarg = sys.argv[2]
    
    int(xarg)
    int(yarg)
except:
    print('Usage: python plot.py <x> <y>')
    print('\t<x> and <y> are the x and y coordinates of the data to be plotted')

print(xarg)
print(yarg)

# Open and parse file.
with open('2D-01-3-samples-per-point.log', 'r') as f:

    line = f.readline();
    while line:
        if 'Ready' in line:
            
            xy = []
            
            # Get coordinates.
            line = re.sub('[^0-9]', ' ', line)
            for i in line.split():
                xy.append(i)
            x = xy[0]
            y = xy[1]

            if (x != xarg or y != yarg):
                continue
           
            # Skip fluff line
            line = f.readline()
            
            # Line 3: Data
            line = f.readline()
            if 'init:' in line:
                
                # Remove non-digit characters.
                line = re.sub('[^0-9]', ' ', line)

                # List to hold data for each mic.
                m0data = []
                m1data = []
                m2data = []
                m3data = []
                # List to hold string data.
                stringdata = line.split()

                # Remove first (fluff) item from string list.
                del stringdata[0]

                # Generate data list for each mic
                for i in range(0, len(stringdata), 2):
                    value = int(stringdata[i])
                    iterations = int(stringdata[i+1])
                    
                    m0 = (value >> 0) & 1
                    m1 = (value >> 1) & 1
                    m2 = (value >> 2) & 1
                    m3 = (value >> 3) & 1
                    
                    for i in range(0, iterations):
                        m0data.append(m0)
                        m1data.append(m1)
                        m2data.append(m2)
                        m3data.append(m3)

                # List to hold x-values
                xvals = []
                for i in range(0, 1024):
                    xvals.append(i)
            
                # Create plot
                fig, axs = plt.subplots(4)
                axs[0].plot(xvals, m0data)
                axs[0].set_title('Mic 0')
                axs[1].plot(xvals, m1data)
                axs[1].set_title('Mic 1')
                axs[2].plot(xvals, m2data)
                axs[2].set_title('Mic 2')
                axs[3].plot(xvals, m3data)
                axs[3].set_title('Mic 3')
                fig.suptitle('Microphone Data for Point (' + x + ',' + y + ')')
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                plt.savefig('plots/(' + x + ',' + y + ').png')
                
                break

        line = f.readline()