clc,clear,close all

name = '18830';
latitude = -22.132705; % degrees north (- = south)
longitude = 133.413717; % degrees east (- = west)

data = load('titree1.txt'); % gage
k=30:length(data);
t1 = data(k,1) + 693960; % excel to matlab date
w1 = data(k,2);
T1 = data(k,3);
clear data

w1 = fixdata(w1,0.3); % Remove steps/spikes > 0.3 cm

data = load('titree2.txt'); % absolute & barometric
k=20:length(data);
t2 = data(k,1) + 693960;
T2 = data(k,2);
h2 = data(k,3);
bp = data(k,4);
clear data

bp = fixdata(bp,1);
w2 = fixdata(h2-bp,1); % convert to gage pressure
bp = bp - 1033.227; % subtract mean air pressure

tide = ['O_1   ';'P_1   ';'S_1   ';'K_1   ';'\phi_1';'\psi_1';'N_2   ';'M_2   ';'S_2   ';'K_2   '];
period = [25.819,24.066,24,23.934,23.869,23.804,12.658,12.421,12,11.967];
pj=1:10; 
p=period(pj); %'/24; tide=tide(pj,:);
np = length(p); 
omega = 2*pi./p; 

dlmwrite('Matlab_files/period.arr',period,'delimiter','','precision','%14.6e\n','newline', 'pc');
dlmwrite('Matlab_files/p.arr',p,'delimiter','','precision','%14.6e\n','newline', 'pc');
dlmwrite('Matlab_files/omega.arr',omega,'delimiter','','precision','%14.6e\n','newline', 'pc');

dt = 1/24;
t = (t2(1):dt:t2(length(t2)))';
y = interp1(t1,w1,t);
x = interp1(t2,bp,t);
dx = diff(x)/dt;
dy = diff(y)/dt;

dlmwrite('Matlab_files/t.arr',t,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/x.arr',x,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/y.arr',y,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/dx.arr',dx,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/dy.arr',dy,'delimiter','','precision','%14.6e','newline', 'pc');

n = length(dx); 
nn=(1:n)'; 
lag = [0:24]; % was 240
nm = length(lag);
v = zeros(n,nm);
for i=1:nm;
	j = lag(i); 
    k = 1:(n-j);
	v(j+k,i) = dx(k);
end 

u1 = zeros(n,np); u2=u1; 
for i=1:np
	tau = omega(i)*t(nn);
	u1(:,i) = cos(tau); 
    u2(:,i) = sin(tau);
end

X = [v,u1,u2]; Z = [ones(n,1),X];

dlmwrite('Matlab_files/tau.arr',tau,'delimiter',' ','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/v.mat',v,'delimiter',' ','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/u1.mat',u1,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/u2.mat',u2,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/X.mat',X,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/Z.mat',Z,'delimiter','','precision','%14.6e','newline', 'pc');

c = Z\dy; nc = length(c); 

dlmwrite('Matlab_files/c.arr',c,'delimiter','','precision','%14.6e','newline', 'pc');

py = y - dt*cumsum([0;X*c(2:nc)]);
oerror = std(dy); 
perror = std(dy - Z*c);

dlmwrite('Matlab_files/py.arr',py,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/oerror.arr',oerror,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/perror.arr',perror,'delimiter','','precision','%14.6e','newline', 'pc');

trend = c(1);
brf = c(1 + (1:nm)); 
crf = abs(cumsum(brf)); 
k = (nm+1)+(1:np); 
trf = complex(c(k),c(np+k));
mag=abs(trf); phase=angle(trf);

dlmwrite('Matlab_files/trend.arr',trend,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/brf.arr',brf,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/crf.arr',crf,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/trf.arr',trf,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/mag.arr',mag,'delimiter','','precision','%14.6e','newline', 'pc');
dlmwrite('Matlab_files/phase.arr',phase,'delimiter','','precision','%14.6e','newline', 'pc');

figure
lag(1)=.25;
semilogx(lag/24,crf,'.-.'),grid
xlabel('Lag (days)'),ylabel('Cumulative Response Function')
title(['Ti Tree Bore ',name])
%print -dtiff Matlab_files/titree1
saveas(gcf,'Matlab_files/titree1.png')

figure
plot(mag,phase,'r.'),grid
xlabel('Tidal Component Amplitude (cm)')
ylabel('Tidal Component Phase (radians)')
title(['Ti Tree Bore ',name])
for i=1:np, text(mag(i)+.02,phase(i),tide(i,:)); end
%print -dtiff Matlab_files/titree2
saveas(gcf,'Matlab_files/titree2.png')

figure
plot(t,y,t,py-1)
datetick('x','mmmyyyy'),ylabel('Relative Water Level')
legend('Observed','Deconvoluted (offset = -1)')
title(['Ti Tree Bore ',name])
to = datenum(2013,5,1);	yo = 145; 
text(to,yo,['Original standard error = ',num2str(oerror),' cm'])
text(to,yo-1,['Residual standard error = ',num2str(perror),' cm'])
text(to,yo-2,['Error ratio = ',num2str(perror/oerror),' cm/cm'])
text(to,yo-3,['Trend = ',num2str(trend),' cm/day'])
xlim([datenum(2013,4,5),datenum(2013,9,10)]);
%print -dtiff Matlab_files/titree3
saveas(gcf,'Matlab_files/titree3.png')

function y = fixdata(x,jump)

% remove data jumps (steps/spikes)

delta = [0;diff(x)]; 
delta(abs(delta)>jump) = 0; 
y = x(1)+cumsum(delta);
end
