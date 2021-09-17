function ierr = add_density_hycom(f)
% add the density field to hycom glb output
% requires the Matlab seawater toolbox (gsw)

g = 9.807; %gravitational acceleration (needed for pressure computation)

% use gsw_rho ( S (g/kg), T (C), p (dbar)) TEOS10
addpath(genpath('/opt/matlab/gsw'));


% read in salinity, temperature, depth data from HYCOM output
time = ncread(f,'time');
ntimes = numel(time);
depth = ncread(f,'depth');
salt = ncread(f,'salinity');
temp = ncread(f,'water_temp');
ndim = numel(size(temp));
if(ndim ==3)
[nx,ny,nd] = size(temp);
else
    [nx,ny,nd,~] = size(temp);
end



density = zeros(size(salt)); %initialize 3D array to hold in-situ density

% loop over horizontal coordinates and compute density in 1 water column at a time
for n=1:ntimes
    fprintf('density processing time %d\n',n);
    for i=1:nx
        
        for j=1:ny
            
            % s,t,rho from water column from surface to bottom
            if(ndim==3)
                s = squeeze(salt(i,j,:));
                t = squeeze(temp(i,j,:));
            else
                s = squeeze(salt(i,j,:,n));
                t = squeeze(temp(i,j,:,n));
            end
            rho = zeros(size(s));

            rho(1) = gsw_rho(s(1),t(1),0); %surf density
            p = 0; %set pressure to zero at the surface
            
            % proceed down water column, accumulate pressure and estimate rho
            for d=2:nd
                p = p + rho(d-1)*g*(depth(d)-depth(d-1));
                rho(d) = gsw_rho(s(d),t(d),p/10000.);
            end
            
            % set density from water column into 3D density array
            if(ndim==3)
                density(i,j,:) = rho;
            else
                density(i,j,:,n) = rho;
            end
        end
    end
end
% plot one column
% nxc = ceil(nx/2); nyc = ceil(ny/2);
% nxc = 20; nyc = 20;
% subplot(1,3,1);
% plot(squeeze(salt(nxc,nyc,:)),-depth);
% subplot(1,3,2);
% plot(squeeze(temp(nxc,nyc,:)),-depth);
% subplot(1,3,3);
% plot(squeeze(density(nxc,nyc,:)),-depth); hold on;

% dump density to netcdf file - save as a short
% with this scaling/offset should be able to save density form 1020-1060 kg/m^3 to precision of
% 0.0005
nccreate(f,'rho','Dimensions',{'lon' nx , 'lat' ny, 'depth' nd, 'time' ntimes},'Datatype','int16','Format','classic');
ncwriteatt(f,'rho','long_name','in-situ density');
ncwriteatt(f,'rho','units','kgm^-3');
ncwriteatt(f,'rho','_FillValue',-30000);
ncwriteatt(f,'rho','scale_factor',0.001);
ncwriteatt(f,'rho','add_offset',1020.);
ncwriteatt(f,'rho','calc_method','TEOS10-gsw_rho-matlab');
ncwrite(f,'rho',density);

% dump matlab dnum to netcdf file
time = ncread(f,'time');
d1 = datenum(2000,01,01);
dnum = datenum(d1 + time/24.);
nccreate(f,'dnum','Dimensions',{'time' ntimes},'Datatype','single','Format','classic');
ncwriteatt(f,'dnum','long_name','Matlab datenumber');
ncwriteatt(f,'dnum','units','days');
ncwrite(f,'dnum',dnum);

% rho_read = ncread(f,'rho');
% plot(squeeze(density(20,20,:)),-depth,'r'); hold on; plot(squeeze(rho_read(20,20,:)),-depth,'bo');
%
