% download hycom data using daterange and bounding box
% add density field and matlab datenum using add_density_hycom
% optionally manipulate velocity/density fields to create idealized scenarios
% G. Cowles 2021

% date range to download
d1 = datenum(2016,1,1);
d2 = datenum(2016,1,4);
dskip = 1;

% system-dependent wget command
wget = '/opt/local/bin/wget -v ';
ncks = '/opt/local/bin/ncks ';

% hycom experiment location online
%locale = ' ''http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.1?';
locale = ' ''http://ncss.hycom.org/thredds/ncss/GLBv0.08/expt_56.3?';

% variables to include
vars = 'var=surf_el,water_u,water_v,water_temp,salinity';

% output prefix
prefix = 'GOM';

% output directory location
%outdir = '/Volumes/FastDrive/ONR_MUST_LatMix_DATA/OceanParcelsSet';
outdir = '../hycom_data';

% lon/lat bounding box
%bbox = '&north=25&south=5&west=-120&east=-100';
bbox = '&north=43.5&south=43&west=-70&east=-69';

% hour string 
%hrs = ['00';'03';'06';'09';'12';'15';'18';'21'];  %uncomment to have every frame in a different
%file
hrs = ['00']; % uncomment for daily files (8 frames in each for 3-hour HYCOM)
[nhrs,~] = size(hrs);

%time = '&time=2019-01-01T00%3A00%3A00Z';
extra = '&addLatLon=true&accept=netcdf'' ';


for d=d1:dskip:d2  % loop over days
    for h=1:nhrs  % loop over hours in the day
        
        % set the wget command string
        if(nhrs > 1)
            time  = ['&time=' datestr(d,'yyyy-mm-dd') 'T' hrs(h,:) '%3A00%3A00Z' ];
        else
            time_start = ['&time_start=' datestr(d,'yyyy-mm-dd') 'T00'  '%3A00%3A00Z'];
            time_end   = ['&time_end=' datestr(d,'yyyy-mm-dd') 'T21' '%3A00%3A00Z'];
            time = [time_start time_end];
        end

        %%time  = ['&time_start=2010-01-04T00:00:00Z&time_end=2010-01-04T00:00:00Z];
        % /opt/local/bin/wget -v  'http://ncss.hycom.org/thredds/ncss/GLBv0.08/expt_56.3?var=surf_el,water_u,water_v,water_temp,salinity&north=43.5&south=43&west=-70&east=-69&time_start=2016-01-01T00%3A00%3A00Z&time_end=2016-01-03T00%3A00%3A00&addLatLon=true&accept=netcdf' -O testday.nc
        %ofile = [datestr(d,'yyyymmdd') hrs(h,:) '0000' ];
        output = [outdir '/' prefix '_'   datestr(d,'yyyymmdd') hrs(h,:) '0000.nc'];
        if(exist(output,'file')); delete(output); end
        cmd = [ wget  locale vars bbox time extra ' -O '  output];
      
        % execute wget
        system(cmd);
        ctmp = cmd;
        
        % check time
        time = ncread(output,'time');
        fprintf('THIS IS TIME:  %f\n',time);
       
        % change time to unlimited / record dimension and dump to tmp.nc
        if(exist('tmp.nc','file')); delete('tmp.nc'); end
        cmd = [ncks '--mk_rec_dmn time ' output ' tmp.nc']; 
        system(cmd);
        
        % overwrite output file with tmp.nc
        cmd = ['mv tmp.nc ' output];
        system(cmd);
        
     
        
        % add the in situ density and matlab dnum to the hycom file
        add_density_hycom(output);
        
        %%%%%%%%%%%%%%% make fake data for checking %%%%%%%%%%%%%%%%%%%%%%
        filetime = (ncread(output,'time')-140256)/24.;
        lon = ncread(output,'lon');
        lat = ncread(output,'lat');
        u = ncread(output,'water_u');
        depth = ncread(output,'depth');
        % velocity field is crude inertial motion
        if(nhrs > 1)
        tmp = 0.05*ones(size(u)) + .1*cos(2*pi*filetime);
        ncwrite(output,'water_u',tmp);
        tmp = 0.05*ones(size(u)) + .1*sin(2*pi*filetime);
        ncwrite(output,'water_v',tmp);
        fprintf('time %f  v %f \n',filetime,tmp(1,1,1));
        else
        tmp = zeros(size(u));
        tmp2 = zeros(size(u));
        for ii=1:numel(filetime)
            tmp(:,:,:,ii) = 0.05 + .1*cos(2*pi*filetime(ii));
            tmp2(:,:,:,ii) = 0.05 + .1*sin(2*pi*filetime(ii));
        end
        ncwrite(output,'water_u',tmp);
        ncwrite(output,'water_v',tmp2);
        fprintf('time %f  v %f \n',filetime,tmp(1,1,1,:));
        
        % create a linearly continuously stratified density field
        for ii=1:numel(depth)
            fac = (1-depth(ii)/max(depth));
           tmp(:,:,ii,:) = fac*1025 + (1-fac)*1035;
           ncwrite(output,'rho',tmp);
        end
        end
        
        %%%%%%%%%%%%%% END MAKE FAKE DATA %%%%%%%%%%%%%%%%%%%%%%%%
        
    end
end