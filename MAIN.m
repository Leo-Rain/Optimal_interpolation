%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS IS THE MAIN PROGRAM FOR OPTIMAL INTERPOLATION (several measurement points + LAGS)  
% LAST MODIFICATION: 11.09.2017
% Whole ARCTIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% BASIC CONSTANTS
%%%%%%%%%%%%%%%%%%%%%
Numb_obs_file = 10; % Number of the file that will be a validation one <==== !!!
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
beta_u = 4;          % - Ridge regression parameter
beta_v = 4;          % - Ridge regression parameter
%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%
%  Month_of_interest = 9; % Variable parameter 
%%%%%%%%%%%%%%%%%%%%%%%%

for Month_of_interest = 9:9
 
% Month_of_interest = 9

% Field plot parameters
%-------------------------------------------------------------------------
t = 350;           % Time instance when we want to plot a field

%-------------------------------------------------------------------------

files_measurements = dir('Data\Obs\rp5_obs_csv26\*.csv'); 
M  = length(files_measurements) - 1;            % Number of files (points) for observations - one validation point
clear files_measurements;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                     START CALCULATION FOR A MONTH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract WRF results for the month of interest:

% ncdisp('Data/WRF_2014/wrfout201409.nc')
T_sample_wrf = 1;
[WRF_Month] = wrf_month(T_sample_wrf, Start_coord, Month_of_interest); % Functon that processes WRF modelled data

u_WRF_Month = WRF_Month{1};
v_WRF_Month = WRF_Month{2};
time_WRF_Month = WRF_Month{3};

%% Extracting time and time length of the Month time seria:

Month_size = size(time_WRF_Month);  % Obtain dimensions of our data
T_Month = Month_size(1);            % Time length, dt = 1h in the Month
clear Month_size;

%% !!! FOR THE ERROR RELATIVE VARIANCE:

% u_YEAR_wrf_var = var(u_YEAR_wrf,0,3);
% v_YEAR_wrf_var = var(v_YEAR_wrf,0,3);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set coordinates of the observation points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Yn_coord, Xj_indices] = RP5_Coord_and_Indexes(M_all, long,lat);

% Xj_indices_short = Xj_indices - ( (Start_coord - 1) * ones(M_all,2) );

Xj_indices_short = Xj_indices;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRF Month at RP5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xj_u_Month_wrf = corrsp_field_serie(u_WRF_Month, M_all, length(time_WRF_Month), Xj_indices_short); 
Xj_v_Month_wrf = corrsp_field_serie(v_WRF_Month, M_all, length(time_WRF_Month), Xj_indices_short);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           MEASUREMENTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract measured RP5 data for U and V components, and upsample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
Yn_u_Month = Observed_data_processing_u_all (time_WRF_Month, M_all, length(time_WRF_Month));
Yn_v_Month = Observed_data_processing_v_all (time_WRF_Month, M_all, length(time_WRF_Month));

%%
% CELLS -> ARRAYS (measured RP5 data)
Yn_u_Month_c = zeros(M_all,length(time_WRF_Month));
for i = 1:M_all
	Yn_u_Month_c(i,:) = Yn_u_Month{i,4};
end

Yn_v_Month_c = zeros(M_all,length(time_WRF_Month));
for i = 1:M_all
	Yn_v_Month_c(i,:) = Yn_v_Month{i,4};
end

% clear Yn_u_Month
% clear Yn_u_Month

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                        MAIN CALCULATIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Rij and Rjj modification - exclude validation point - CHECK IT !!!

Rij_i = Rij;
size_Rij = size(Rij);
half = size_Rij(3)/2;
Rij_i(:,:,[Numb_obs_file, Numb_obs_file + half ]) = [];

Rjj_i = Rjj;
Rjj_i(:,[Numb_obs_file, Numb_obs_file + half]) = [];
Rjj_i([Numb_obs_file, Numb_obs_file + half],:) = [];

%% KALMAN WEIGHTS

% Exclude verification point from RP5 WRF
 Xj_u_Month_wrf_i =  Xj_u_Month_wrf;
 Xj_u_Month_wrf_i(Numb_obs_file,:) = [];
 
 Xj_v_Month_wrf_i =  Xj_v_Month_wrf;
 Xj_v_Month_wrf_i(Numb_obs_file,:) = [];
  
% Exclude verification point from RP5 measurements
Yn_v_Month_c_i = Yn_v_Month_c;
Yn_v_Month_c_i(Numb_obs_file,:) = [];

Yn_u_Month_c_i = Yn_u_Month_c;
Yn_u_Month_c_i(Numb_obs_file,:) = [];

%%
% Kalman weights estimation - solving the linear system here!
a = ridge(Rjj_i, Rij_i, Xj_u_Month_wrf_i, Xj_v_Month_wrf_i, Yn_u_Month_c_i, Yn_v_Month_c_i, beta_u, beta_v);

%% INNOVATION FIELD

[d_u, d_v] = calc_correction_2(a, Xj_u_Month_wrf_i, Xj_v_Month_wrf_i, Yn_u_Month_c_i, Yn_v_Month_c_i, X, Y, M, T, T_Month); % Main calculation

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Extract RP5 data for the Reference point for U and V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Measured RP5 data at the ref point ( = Numb_obs_file) 
[Yn_u_ref, Ref_point_coord] = Measured_data_processing_u_ref(time_WRF_Month,length(time_WRF_Month), Numb_obs_file);
 Yn_v_ref                   = Measured_data_processing_v_ref(time_WRF_Month,length(time_WRF_Month), Numb_obs_file);

 Ref_point_pos = [0 0]; % initialize x and y of the reference point where we take out time plots 
[Ref_point_pos(1,1),Ref_point_pos(1,2)] = near2(lat,long,Ref_point_coord(1,1),Ref_point_coord(1,2)); % Indices of the ref point

% Ref_point_pos_i = Ref_point_pos - ( (Start_coord - 1) * ones(1,2) );
Ref_point_pos_i = Ref_point_pos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                        RELATIVE ANALYSIS ERROR VARIANCE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% 
% u_YEAR_wrf_var = var(u_YEAR_wrf,0,3);
% v_YEAR_wrf_var = var(v_YEAR_wrf,0,3);
% 
% % Std full
% u_YEAR_wrf_var_full = zeros(Y,X);
% v_YEAR_wrf_var_full = zeros(Y,X);
% 
% Start_coord = 1;
% 
% u_YEAR_wrf_var_full(Start_coord:end,Start_coord:end) = u_YEAR_wrf_var; % U mean for 2014
% v_YEAR_wrf_var_full(Start_coord:end,Start_coord:end) = v_YEAR_wrf_var; % V mean for 2014
% 
% %%
% [P_u,P_v] = Analysis_error(a,Rij_i,X,Y, u_YEAR_wrf_var_full, v_YEAR_wrf_var_full);
% 
% %%
% % PLOT ERROR VARIANCE FOR EACH COMPONENT
% 
% % FOR U:
% figure
% m_map_projections_scalar(P_u,lat,long,lat_lim,lon_lim)
% %plot measurement points:
% for i = 1:M
%     plot_point(Yn_coord(i,1),Yn_coord(i,2),i);
% end
% 
% % % plot reference point: 
% % plot_point_reference(Ref_point_coord(2),Ref_point_coord(1));
% % tune the plot:
% color_max = 1;
% color_min = 0;
% caxis([color_min color_max]);
% colormap(jet)
% % colormap(flipud(colormap))
% c = colorbar('eastoutside');
% title('Relative error variance, U');
% 
% %%
% % FOR V:
% figure
% m_map_projections_scalar(P_v,lat,long,lat_lim,lon_lim)
% %plot measurement points:
% for i = 1:M
%     plot_point(Yn_coord(i,1),Yn_coord(i,2),i);
% end
% 
% % % plot reference point: 
% % plot_point_reference(Ref_point_coord(2),Ref_point_coord(1));
% % tune the plot:
% color_max = 1;
% color_min = 0;
% caxis([color_min color_max]);
% colormap(jet)
% % colormap(flipud(colormap))
% c = colorbar('eastoutside');
% title('Relative error variance, V');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                 FIELDS AND TIME SERIES OF INTEREST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_u_t = d_u(:,:,t);
d_v_t = d_v(:,:,t);

% Correction term addition: WRF_Aug_u(439*479*t) + du(439*479*t)
u_t = u_WRF_Month(:,:,t);
u_t_corrected = u_t + d_u_t;

% Correction term addition: WRF_Aug_v(439*479*t) + dv(439*479*t)
v_t = v_WRF_Month(:,:,t);
v_t_corrected = v_t + d_v_t;

% Now correction term addition for the whole field and whole time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES TO SAVE:
u_corrected = u_WRF_Month + d_u; % corrected field for component u
v_corrected = v_WRF_Month + d_v; % corrected field for component v

%%
% % Save them to *.nc files:
% cd('Data\WRF_corrected'); 
% 
%     nccreate( ['corrected_wind_', num2str(Month_of_interest), '.nc'],'u_10m_corrected',   'Datatype', 'double', 'Format', '64bit', 'Dimensions',{'dimY',Y,'dimX',X,'time',length(time_WRF_Month)} );
%     nccreate( ['corrected_wind_', num2str(Month_of_interest), '.nc'],'v_10m_corrected',   'Datatype', 'double', 'Format', '64bit', 'Dimensions',{'dimY',Y,'dimX',X,'time',length(time_WRF_Month)} );
%     nccreate( ['corrected_wind_', num2str(Month_of_interest), '.nc'],'Time',              'Datatype', 'double', 'Format', '64bit', 'Dimensions',{'time',length(time_WRF_Month), 'dim', 1} );
%     
%     ncwrite(['corrected_wind_', num2str(Month_of_interest), '.nc'],'u_10m_corrected', u_corrected);
%     ncwrite(['corrected_wind_', num2str(Month_of_interest), '.nc'],'v_10m_corrected', v_corrected);
%     ncwrite(['corrected_wind_', num2str(Month_of_interest), '.nc'],'Time', time_WRF_Month);
%     
% % ncdisp('corrected_wind_9.nc')
% cd('C:\Users\agusarov\Desktop\WRF_BA_5\Code'); % Go to dir with files on rp5 of interest

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Magnitude wind calculation
A_t = sqrt((u_t.^2) + (v_t.^2)); % for original field at t
A_corrected_t = sqrt((u_t_corrected.^2) + (v_t_corrected.^2)); % for corrected field at t

% Correction field for the magnitude 
d = sqrt((d_u(:,:,:).^2) + (d_v(:,:,:).^2)); % Whole field, whole time
d_t = d(:,:,t); % Whole field, at t

d_u_ref_t = d_u(Ref_point_pos_i(1), Ref_point_pos_i(2),:); % Correction to U time serie at the ref point
d_u_ref_t = reshape(d_u_ref_t ,1,[]);                      % Make a linear array out of it (just reshape)
%------------------------------
u_ref = u_WRF_Month(Ref_point_pos_i(1), Ref_point_pos_i(2),:); % Time-serie U at the ref. point
    u_ref = reshape( u_ref ,1,[]); % Make a linear array out of it (just reshape)

u_ref_corrected = u_corrected(Ref_point_pos_i(1), Ref_point_pos_i(2),:); % corrected field for component U in ref. point
    u_ref_corrected = reshape( u_ref_corrected ,1,[]); % Make a linear array out of it
%------------------------------
v_ref = v_WRF_Month(Ref_point_pos_i(1), Ref_point_pos_i(2),:); % Time-serie V at the ref. point
    v_ref = reshape( v_ref ,1,[]); % Make a linear array out of it

v_ref_corrected = v_corrected(Ref_point_pos_i(1), Ref_point_pos_i(2),:); % corrected field for component v in ref. point
    v_ref_corrected = reshape( v_ref_corrected ,1,[]); % Make a linear array out of it

end 
    
    

