import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

mode = 1 ## (mode 1: directly plot, mode 2: save to movie)
f = open("Data_ParticlePosition")
lines = f.readlines()
f.close()

Test = np.array(list(map(float,lines[1].split())))
N_num = Test.shape[0]
ParN = int((N_num-1)/3)
T_num = len(lines)-1

RawData = np.zeros((T_num, N_num))

for i in range(T_num):
    RawData[i,:] = np.array(list(map(float,lines[i+1].split())))

PlotData = np.swapaxes(RawData,0,1)

def GetData(input,index):
    gdata = np.zeros((3,T_num))
    for T in range(T_num):
        for D in range(3):
            gdata[D,T] = PlotData[3*index+1+D, T]
    return gdata

data = [GetData(PlotData, index) for index in range(ParN)]


pvs = 0.01  #particle size of view
L = 1

def update_lines(num, dataLines, lines):
    for line, data in zip(lines, dataLines):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
    return lines

# Attaching 3D axis to the figure

# Set up formatting for the movie files

fig = plt.figure()
ax = p3.Axes3D(fig)

# Fifty lines of random 3-D lines

# Creating fifty line objects.
# NOTE: Can't pass empty arrays into 3d version of plot()
lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

# Setting the axes properties
ax.set_xlim3d([0.0, L])
ax.set_xlabel('X')

ax.set_ylim3d([0.0, L])
ax.set_ylabel('Y')

ax.set_zlim3d([0.0, L])
ax.set_zlabel('Z')

ax.set_title('Particle Trajectory')

# Creating the Animation object

line_ani = animation.FuncAnimation(fig, update_lines, T_num, fargs=(data, lines), interval=5, blit=False)

if mode == 1:
    plt.show()
else:
    line_ani.save('mymovie.mp4') ##https://morvanzhou.github.io/tutorials/data-manipulation/plt/5-1-animation/
