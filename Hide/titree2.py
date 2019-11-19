
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
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
tide = ['O$_1$   ', 'P$_1$   ', 'S$_1$   ', 'K$_1$   ', '$\phi_1$',
        '$\psi_1$', 'N$_2$   ', 'M$_2$   ', 'S$_2$   ', 'K$_2$   ']
np.savetxt('Data/tide.dat', tide, fmt='%s')
period = [  25.819,   24.066,   24.000,   23.934,   23.869, 
            23.804,   12.658,   12.421,   12.000,   11.967]
p     = np.array(period) 
np.savetxt('Data/p.dat', p)
NP    = len(p)
omega = 2.*np.pi/p # angular frequency calculated for specified periods
np.savetxt('Data/omega.dat', omega)

dt = 1./24.

data = np.loadtxt('Data/titree1.txt', delimiter='\t')
#data = np.loadtxt('Data/'+name+'.txt', delimiter='\t')
k  = range(30-1, len(data))
t1 = data[k,0]
w1 = data[k,1]          
T1 = data[k,2]          
w1 = fixdata(w1, 0.3)
np.savetxt('Data/t1.dat',  t1)
np.savetxt('Data/w1.dat',  w1)
np.savetxt('Data/t1b.dat', T1)

data = np.loadtxt('Data/titree2.txt', delimiter='\t')
k  = range(20-1, len(data))
t2 = data[k,0]
T2 = data[k,1]
h2 = data[k,2]
bp = data[k,3]
bp = fixdata(bp, 1.)
w2 = fixdata(h2-bp, 1.) # convert to gauge pressure
bp = bp-1033.227 # subtract mean air pressure
np.savetxt('Data/t2.dat',  t2)
np.savetxt('Data/T2b.dat', T2)
np.savetxt('Data/h2.dat',  h2)
np.savetxt('Data/w2.dat',  w2)
np.savetxt('Data/bp.dat',  bp)

t = np.arange(t2[0], t2[-1]+dt, dt)
y = np.interp(t, t1, w1)
x = np.interp(t, t2, bp) 
dx = np.diff(x)/dt
dy = np.diff(y)/dt
np.savetxt('Data/t.dat',  t)
np.savetxt('Data/x.vec',  x)
np.savetxt('Data/y.vec',  y)
np.savetxt('Data/dx.dat', dx)
np.savetxt('Data/dy.dat', dy)

n   = len(dx)
nn  = range(n) 
lag = range(240+1)
np.savetxt('Data/nn.dat',  nn)
np.savetxt('Data/lag.dat', lag)
nm = len(lag)
v = np.zeros([n, nm])
for i in range(nm):
    j = lag[i]
    k = np.arange(n-j)
    v[j+k, i] = dx[k] 
np.savetxt('Data/v.mat', v)

u1 = np.zeros([n, NP])
u2 = u1.copy()
for i in range(NP):
    tau = omega[i]*t[nn]
    u1[:,i] = np.cos(tau)
    u2[:,i] = np.sin(tau)

X = np.hstack([v, u1, u2])
Z = np.hstack([np.ones([n,1]), X])
np.savetxt('Data/X.mat', X)
np.savetxt('Data/Z.mat', Z)

c  = np.linalg.lstsq(Z, dy, rcond=None)[0]
nc = len(c)
np.savetxt('Data/c.dat', c)

py = y-dt*np.concatenate([[0.], np.cumsum(np.dot(X, c[1:nc]))])
oerror = np.std(dy)
perror = np.std(dy-np.dot(Z,c))
np.savetxt('Data/py.dat', py)

trend = c[0]
brf   = c[0+np.arange(nm)]
crf   = np.abs(np.cumsum(brf))
k     = np.concatenate([np.atleast_1d([nm+1]), np.arange(NP)])
trf   = [a+1j*b for a,b in zip(c[k], c[NP+k])]
mag   = np.abs(trf)
phase = np.angle(trf)


t = [datetime.utcfromtimestamp(i) for i in (t-25569.)*86400.]

f,s = plt.subplots()
lag[1] = 0.25
s.set_title('Ti Tree Bore '+name, fontsize=10, fontweight='bold')
plt.semilogx(np.array(lag)/24., crf, 'ko-', ms=4, mec='k', mfc='none', mew=1.)
s.set_xlabel('Lag (days)')
s.set_ylabel('Cumulative response function')
#s.set_xlim(1e-2, 1)
#s.set_ylim(0.00, 0.07)
for i in ['top', 'right']:
    s.spines[i].set_visible(False)
s.grid(which='major', axis='both', c=(194./255., 194./255., 194./255.), ls='-', 
       lw=0.5)                            
s.grid(which='minor', axis='x', c=(194./255., 194./255., 194./255.), ls='-', 
       lw=0.5)                            
plt.tight_layout()
plt.savefig('Figs/titree1.png', dpi=300)
plt.close(f)

f,s = plt.subplots()
s.set_title('Ti Tree Bore '+name, fontsize=10, fontweight='bold')
s.plot(mag, phase, 'ro')
s.set_xlabel('Tidal component amplitude (cm)')
s.set_ylabel('Tidal component phase (radians)')
for i in range(NP):
    s.text(mag[i]+5e-4, phase[i]+5e-2, tide[i])
#s.set_xlim(0.00, 0.05)
#s.set_ylim(-4., 3.)
for i in ['top', 'right']:
    s.spines[i].set_visible(False)
s.grid(which='major', axis='both', c=(194./255., 194./255., 194./255.), ls='-', 
       lw=0.5)                            
plt.tight_layout()
plt.savefig('Figs/titree2.png', dpi=300)
plt.close(f)

f,s = plt.subplots()
s.set_title('Ti Tree Bore '+name, fontsize=10, fontweight='bold')
s.plot(t, y,     'k-', label='Observed')
s.plot(t, py-1., 'r-', label='Deconvoluted (offset = -1)')
s.set_ylabel('Relative water level (cm)')
s.legend(fancybox=False)
to = datetime(2013, 5, 1)
yo = 145.25 
s.text(to, yo,   'Original standard error = '+str('%.4f'% oerror)+' cm')
s.text(to, yo-1, 'Residual standard error = '+str('%.4f'% perror)+' cm')
s.text(to, yo-2, 'Error ratio = '+str('%.4f'% (perror/oerror))+' cm/cm')
s.text(to, yo-3, 'Trend = '+str('%.4f'% trend)+' cm/day')
s.xaxis.set_major_formatter(mpl.dates.DateFormatter('%b\n%Y'))
#s.set_xlim(datetime(2013, 4, 5), datetime(2013, 9, 10))
#s.set_ylim(140., 158.)
for i in ['top', 'right']:
    s.spines[i].set_visible(False)
s.grid(which='major', axis='both', c=(194./255., 194./255., 194./255.), ls='-', 
       lw=0.5)                            
plt.tight_layout()
plt.savefig('Figs/titree3.png', dpi=300)
plt.close(f)
