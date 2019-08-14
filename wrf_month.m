function[wrf_months] = wrf_month(T_sample, Start_coord, Month_of_interest)

N_months  = 1;                         % Number of months 
files_wrf = dir('Data\WRF_2014\*.nc'); 
cd('Data\WRF_2014');                    % Go to dir with files on rp5 of interest

wrf_months = cell(N_months,3);          % N_months x (u,v,xtime)


    file = files_wrf(Month_of_interest,:);

    wrf_months{1,1} = ncread(file.name,'u_10m', [Start_coord Start_coord 1],[Inf Inf Inf], [1 1 T_sample]);              
    wrf_months{1,2} = ncread(file.name,'v_10m', [Start_coord Start_coord 1],[Inf Inf Inf], [1 1 T_sample]);      
    wrf_months{1,3} = ncread(file.name,'XTIME'); % Extract time dimension from WRF


    time_conv = wrf_months{1,3}*60 + 6*60*60; 
%     time_conv = wrf_months{1,3}*60; 
    time_conv = datetime(time_conv,'ConvertFrom','epochtime','Epoch','2012-01-01');
    
cd('C:\Users\agusarov\Desktop\WRF_BA_5\Code'); % Go back to home dir

