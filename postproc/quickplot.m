clear all; close all;
f = 'isopycnal_float_multiframe.nc';
z = ncread(f,'z');
rho = ncread(f,'rho');
time = ncread(f,'time');
w = ncread(f,'w');
tmp1 = ncread(f,'tmp1');
lon = ncread(f,'lon');
lat = ncread(f,'lat');

plot(time(:,1)/86400.,rho(:,1),'r'); hold on; plot(time(:,2)/86400.,rho(:,2),'b');