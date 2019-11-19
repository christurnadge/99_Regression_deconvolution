
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from datetime import datetime
from fixdata import fixdata

mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams[ 'font.sans-serif'  ] = 'Calibri'
mpl.rcParams[ 'font.size'        ] = 10
mpl.rcParams[ 'mathtext.default' ] = 'regular'          
mpl.rcParams[ 'xtick.direction'  ] = 'out'
mpl.rcParams[ 'ytick.direction'  ] = 'out'       
mpl.rcParams[ 'lines.linewidth'  ] = 1.0     

name = str(int(np.loadtxt('Data/name.txt')))
latitude  = -22.132705 # degrees north (- = south)
longitude = 133.413717 # degrees east  (- = west)

data = np.loadtxt('Data/titree1.txt', delimiter='\t')
#data = np.loadtxt('Data/'+name+'.txt', delimiter='\t')
k  = range(30-1, len(data))
t1 = (np.array(data[k,0])-25569.)*86400.
w1 = data[k,1] 
T1 = data[k,2]
w1 = fixdata(w1, 0.3)
np.savetxt('Data/t1.dat',  t1)
np.savetxt('Data/w1.dat',  w1)
np.savetxt('Data/T1b.dat', T1)

data = np.loadtxt('Data/titree2.txt', delimiter='\t') 
k  = range(20-1, len(data))
t2 = (np.array(data[k,0])-25569.)*86400.
T2 = data[k,1]
h2 = data[k,2]
bp = data[k,3]
bp = fixdata(bp, 1.)
w2 = fixdata(h2-bp, 1.) # convert to gauge pressure
bp = bp-1033.227 # subtract mean air pressure

data = np.loadtxt('Data/titree3.txt', delimiter='\t')
k  = range(len(data))
t3 = [datetime.utcfromtimestamp(i) for i in (np.array(data[k,0])-25569.)
                                                                *86400.]
w3 = 100.*data[k,2] # convert m to cm
np.savetxt('Data/t2.dat',  t2)
np.savetxt('Data/T2b.dat', T2)
np.savetxt('Data/h2.dat',  h2)
np.savetxt('Data/w2.dat',  w2)
np.savetxt('Data/bp.dat',  bp)
np.savetxt('Data/w3.dat',  w3)

t1 = [datetime.utcfromtimestamp(i) for i in t1]
t2 = [datetime.utcfromtimestamp(i) for i in t2]
      
f,s = plt.subplots()
s.set_title('Ti Tree Bore '+name, fontsize=10, fontweight='bold')

s.plot(t1, T1, 'k-', label='Gauge')

s.plot(t2, T2, 'r-', label='Absolute')

s.set_ylabel('Temperature ($^o$C)')
s.legend(loc=2, fancybox=False)
s.xaxis.set_major_formatter(mpl.dates.DateFormatter('%b\n%Y'))
#s.set_xlim(datetime(2012, 6, 15), datetime(2013, 9, 15))
#s.set_ylim(24.5, 29.0)
for i in ['top', 'right']:
    s.spines[i].set_visible(False)
s.grid(which='major', axis='both', c=(194./255., 194./255., 194./255.), ls='-', 
       lw=1.0)                            
plt.tight_layout()
plt.savefig('Figs/temp.png', dpi=300)
plt.close(f)


f,s = plt.subplots()
s.set_title('Ti Tree Bore '+name, fontsize=10, fontweight='bold')
s.plot(t1, w1-190., 'k-', label='Vented')
#s.set_xlim(datetime(2012, 6, 15), datetime(2013, 9, 15))
s.plot(t2, w2-105., 'r-', label='Gauge')
s.plot(t3, w3-75.,  'b-', label='Capacitance')
s.set_ylabel('Relative water level (cm)')

s.legend([Line2D([0],[0], color='k', lw=1.0),
          Line2D([0],[0], color='r', lw=1.0),
          Line2D([0],[0], color='b', lw=1.0)],
         ['Vented', 'Gauge', 'Capacitance'], loc=1, fancybox=False)

for i in ['top', 'right']:
    s.spines[i].set_visible(False)
s.grid(which='major', axis='both', c=(194./255., 194./255., 194./255.), ls='-', 
       lw=1.0)                            
s.xaxis.set_major_formatter(mpl.dates.DateFormatter('%b\n%Y'))
#s.set_ylim(-50., 0.05)
plt.tight_layout()
plt.savefig('Figs/level.png', dpi=300)
plt.close(f)
