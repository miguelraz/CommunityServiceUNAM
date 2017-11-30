clear all;

N=256;
NEST=20;     %11;
escala=10^2*1; % 4;
offset=49/19;  % 60/19;  % 49/19; %2/19; %(5-0.25)/19;  %164.8;
offini=1;  %-30; 
h=25;  ts=10-0.39;   xini=0; xfin=50; alfa=2.5; beta=1.25;
figure; 

%% Vertical Force
u=load('SINTSxz.RES');
for j=1:NEST;
    ti=N*(j-1)+1;
    tf=N*j;
    plot(u(ti:tf,1)+100*0, u(ti:tf,2)*escala+offset*(j-1)+offini,'b','linewidth',1.5);
    grid on;
    hold on;
end;
v=load('SINTSyz.RES');
for j=1:NEST;
    ti=N*(j-1)+1;
    tf=N*j;
    plot(v(ti:tf,1)+100*1, v(ti:tf,2)*escala+offset*(j-1)+offini,'k','linewidth',1.5);
    grid on;
    hold on;
end;
w=load('SINTSzz.RES');
for j=1:NEST;
    ti=N*(j-1)+1;
    tf=N*j;
    plot(w(ti:tf,1)+100*2, w(ti:tf,2)*escala+offset*(j-1)+offini,'r','linewidth',1.5);
    grid on;
    hold on;
end;

x=xini:0.1:xfin; 
taop=ts+(1/alfa)*sqrt(h^2+x.^2);
plot(taop, x, '--b', 'linewidth', 1);
hold on;
taos=ts+(1/beta)*sqrt(h^2+x.^2);
plot(taos, x, '--r', 'linewidth', 1);
hold on;
plot(taop+100, x, '--b', 'linewidth', 1);
hold on;
plot(taos+100, x, '--r', 'linewidth', 1);
hold on;
plot(taop+200, x, '--b', 'linewidth', 1);
hold on;
plot(taos+200, x, '--r', 'linewidth', 1);
hold on;

%        x = 0:10;  y = sin(x);
%        xx = 0:.25:10;
%        yy = spline(x,y,xx);
%        plot(x,y,'o',xx,yy)
% for j=1:NEST;
%     ti=N*(j-1)+1;
%     tf=N*j;
%     uux=u(ti:0.25:tf,1); uuy=spline(u(ti:tf,1),u(ti:tf,2),uux);
%     wwx=w(ti:0.25:tf,1); wwy=spline(w(ti:tf,1),w(ti:tf,2),wwx);
%     %plot(uux+250,uuy*escala+offset*(j-1)+offini,'r','linewidth',1.5);
%     %plot(wwx+250,wwy*escala+offset*(j-1)+offini,'r','linewidth',1.5);
%     plot(uuy*500+275, wwy*escala+offset*(j-1)+offini,'k','linewidth',1);
%     grid on;
%     hold on;
% end;


xlabel('Time (sec)', 'fontsize', 14)
ylabel('Offset (km)' , 'fontsize', 14) 
title('Seismograms Giz, i=x,y,z', 'fontsize', 16)
axis([0 300  -5  52]);  % x, y, z
%axis([  0 10  -1.2  1.2]); % x
%axis([10 20  -1.2  1.2]);  % y
%axis([20 30  -1.2  1.2]);  % z

figure; 

%% Horizontal Force  Fy
f=load('SINTSxy.RES');
for j=1:NEST;
    ti=N*(j-1)+1;
    tf=N*j;
    plot(f(ti:tf,1)+100*0, f(ti:tf,2)*escala+offset*(j-1)+offini,'b','linewidth',1.5);
    grid on;
    hold on;
end;
f=load('SINTSyy.RES');
for j=1:NEST;
    ti=N*(j-1)+1;
    tf=N*j;
    plot(f(ti:tf,1)+100*1, f(ti:tf,2)*escala+offset*(j-1)+offini,'k','linewidth',1.5);
    grid on;
    hold on;
end;
f=load('SINTSzy.RES');
for j=1:NEST;
    ti=N*(j-1)+1;
    tf=N*j;
    plot(f(ti:tf,1)+100*2, f(ti:tf,2)*escala+offset*(j-1)+offini,'r','linewidth',1.5);
    grid on;
    hold on;
end;

plot(taop, x, '--b', 'linewidth', 1);
hold on;
plot(taos, x, '--r', 'linewidth', 1);
hold on;
plot(taop+100, x, '--b', 'linewidth', 1);
hold on;
plot(taos+100, x, '--r', 'linewidth', 1);
hold on;
plot(taop+200, x, '--b', 'linewidth', 1);
hold on;
plot(taos+200, x, '--r', 'linewidth', 1);

xlabel('Time (sec)', 'fontsize', 14)
ylabel('Offset (km)' , 'fontsize', 14) 
title('Seismograms Giy, i=x,y,z', 'fontsize', 16)
axis([0 300  -5 52]);  % x, y, z
%axis([  0 10  -1.2  1.2]); % x
%axis([10 20  -1.2  1.2]);  % y
%axis([20 30  -1.2  1.2]);  % z

figure;

%% Horizontal Force   Fx
f=load('SINTSxx.RES');
for j=1:NEST;
    ti=N*(j-1)+1;
    tf=N*j;
    plot(f(ti:tf,1)+100*0, f(ti:tf,2)*escala+offset*(j-1)+offini,'b','linewidth',1.5);
    grid on;
    hold on;
end;
f=load('SINTSyx.RES');
for j=1:NEST;
    ti=N*(j-1)+1;
    tf=N*j;
    plot(f(ti:tf,1)+100*1, f(ti:tf,2)*escala+offset*(j-1)+offini,'k','linewidth',1.5);
    grid on;
    hold on;
end;
f=load('SINTSzx.RES');
for j=1:NEST;
    ti=N*(j-1)+1;
    tf=N*j;
    plot(f(ti:tf,1)+100*2, f(ti:tf,2)*escala+offset*(j-1)+offini,'r','linewidth',1.5);
    grid on;
    hold on;
end;

plot(taop, x, '--b', 'linewidth', 1);
hold on;
plot(taos, x, '--r', 'linewidth', 1);
hold on;
plot(taop+100, x, '--b', 'linewidth', 1);
hold on;
plot(taos+100, x, '--r', 'linewidth', 1);
hold on;
plot(taop+200, x, '--b', 'linewidth', 1);
hold on;
plot(taos+200, x, '--r', 'linewidth', 1);

xlabel('Time (sec)', 'fontsize', 14)
ylabel('Offset (km)' , 'fontsize', 14) 
title('Seismograms Gix, i=x,y,z', 'fontsize', 16)
axis([0 300  -5 52]);  % x, y, z
%axis([  0 10  -1.2  1.2]); % x
%axis([10 20  -1.2  1.2]);  % y
%axis([20 30  -1.2  1.2]);  % z

