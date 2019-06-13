import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

PlotMode = 1 ## (mode 1: directly plot, mode 2: save to movie)
ViewMode = 2 ## (mode 1: with trajactory, mode 2: without trajactory)
ps = 1  #particle size of view, only for ViewMode = 2
L = 14    ## size of the BoxSize
filename = "Data_ParticlePosition"
delay = 0.1  ###delay between the nearby frames

## Collect Data

f = open(filename)
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


fig = plt.figure()
ax = p3.Axes3D(fig)
text  = ax.text2D( 0.08, 0.08, '', fontsize=16, color='black', ha='center', va='center' )
if ViewMode == 1:
    lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]
else:
    lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1],'ro',ms=ps)[0] for dat in data]

def update_lines(num, dataLines, lines ):
    text.set_text( 't = %6.3f d' % PlotData[0, num] )
    for line, data in zip(lines, dataLines):
        # NOTE: there is no .set_data() for 3 dim data...
        if ViewMode == 1:
            line.set_data(data[0:2, :num])
            line.set_3d_properties(data[2, :num])
        else:
            line.set_data(data[0:2, num-1:num])
            line.set_3d_properties(data[2, num-1:num])
    return lines

# Setting the axes properties

ax.set_xlim3d([0.0, L])
#ax.set_xlim3d([15, 25])
ax.set_xlabel('X')
ax.set_ylim3d([0.0, L])
#ax.set_ylim3d([15, 25])
ax.set_ylabel('Y')
ax.set_zlim3d([0.0, L])
#ax.set_zlim3d([15, 25])
ax.set_zlabel('Z')
ax.set_aspect( 'equal' )
ax.set_title('Asteroid')

# Creating the Animation object

line_ani = animation.FuncAnimation(fig, update_lines, T_num, fargs=(data, lines), interval=delay, blit=False, repeat = False)
if PlotMode == 1:
    plt.show()
else:
    line_ani.save('PM_movie.mp4') ##https://morvanzhou.github.io/tutorials/data-manipulation/plt/5-1-animation/
