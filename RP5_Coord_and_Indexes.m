function [Yn_coord, Xj_indices] = RP5_Coord_and_Indexes(M, long,lat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set coordinates of the observation points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Yn_coord   = zeros(M,2);   % 1st column = ALT,     2nd column = LONG
Xj_indices = zeros(M,2);   % 1st column = X index, 2nd column = Y index

files_measurements = dir('Data\Obs\rp5_obs_csv26\*.csv'); 
cd('Data\Obs\rp5_obs_csv26'); % Go to dir with files on rp5 of interest
% files_measurements(Numb_obs_file,:)=[];

% Extracting coordinates of RP5
for g = 1:M
    
    file = files_measurements(g,:);
    fileID = fopen(file.name,'r');
           
    C = textscan(fileID,'%s %s %f %f %f %f %f %f %f %f %f %f %s %f %f',...
                        'Delimiter',',','EmptyValue',0,'HeaderLines',1);
    lat_rp_5 = C{1,10};
    lon_rp_5 = C{1,11};
    
    Yn_coord(g,2) = lat_rp_5(1);
    Yn_coord(g,1) = lon_rp_5(1);
    
    fclose(fileID);  
end

 cd('C:\Users\Asus\Desktop\WRF_BA_5'); % Go back to the home dir

% Find indexes according to extracted points' coordinates:
for k = 1:M
    [Xj_indices(k,1), Xj_indices(k,2)] = near2(long,lat,Yn_coord(k,1),Yn_coord(k,2)); 
end
