function [Yn_u, Ref_point_coord] = Measured_data_processing_u_ref (time,N, Numb_obs_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract measured data for U component and upsample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a cell of time series for real measurements:
Yn_u = cell(1,5);

% start and end time in original WRF file
time_start = time(1)  *60 + 6*60*60; % 81475200
time_end   = time(end)*60 + 6*60*60; % 84150000
Subset1 = datetime(time_start,'ConvertFrom','epochtime','Epoch','2012-01-01'); % Start of August {'01-Aug-2014 00:00:00'}
Subset2 = datetime(time_end,  'ConvertFrom','epochtime','Epoch','2012-01-01');   % End of August {'31-Aug-2014 23:00:00'}

% convert 'time' from original WRF (minutes from 2012-01-01) into datetime format:
time_conv = time*60 + 6*60*60; 
time_conv = datetime(time_conv,'ConvertFrom','epochtime','Epoch','2012-01-01');

files = dir('Data\Obs\rp5_obs_csv26\*.csv'); 
files = files(Numb_obs_file,:);
cd('Data\Obs\rp5_obs_csv26')% Go to dir with measurement files

g = 1;
for file = files'
    fileID = fopen(file.name,'r');
        
    C = textscan(fileID,'%{yyyy-MM-dd HH:mm:ss}D %s %f %f %f %f %f %f %f %f %f %f %s %f %f',...
                          'Delimiter',',','EmptyValue',0,'HeaderLines',1);
    U_component = C{1,14};
    Time = C{1,1};
    
    Ref_point_coord = zeros(1,2);
    Ref_point_coord(1,1) = C{1,10}(1);
    Ref_point_coord(1,2) = C{1,11}(1);
    
    
    Superset = datetime(Time,'InputFormat','yyyy-MM-dd HH:mm:ss');
try  
    % Calc index of time start for measured time-serie
    Index1 = datefind(Subset1, Superset);
    if isempty(Index1) == 1 
    i = 1;
    while isempty(Index1) == 1 
        Index1 = datefind(Subset1, Superset,i);
        i = i + 1;
    end
    
    % Index1 = Index1(1);
    [~,ind1] = min( abs(datenum(Subset1)-datenum( Superset(Index1) )) );
    Index1 = Index1(ind1); 
    while Time(Index1) < Subset1 % in case if found date is in less than a target month
        i = 1;
        Index1 = Index1 + i;
        i = i + 1;
    end
    end
    
    % Calc index of time end for measured time-serie
    Index2 = datefind(Subset2, Superset);
    if isempty(Index2) == 1 
    i = 1;
    while isempty(Index2) == 1 
        Index2 = datefind(Subset2, Superset,i);
        i = i + 1;
    end
     % Index2 = Index2(end);
     [~,ind2] = min( abs(datenum(Subset2)-datenum( Superset(Index2) )) );
     Index2 = Index2(ind2); 
     while Time(Index2) > Subset2 % in case if found date is in less than a target month
        i = 1;
        Index2 = Index2 - i;
        i = i + 1;
     end
    end  
    
    U_i = U_component(Index1:Index2);
    T_i = Time(Index1:Index2);
    Yn_u{g,1} = rot90(U_i); % Time serie itself
    Yn_u{g,2} = rot90(T_i); % Time instances for this time-serie 
            
    if isempty(Yn_u{g,1}) == 1 
        Yn_u{g,1} = zeros(1,length(time));
        Yn_u{g,2} = rot90(time_conv);
    end
    
    t1 = datevec(Yn_u{g,2}); % Time atributes of the observed array in timevec format
    for i = 1:length(t1)
        Yn_u{g,3}(i) = etime(t1(i,:),datevec(Subset1))/3600; % elapsed time in hours since start of the month
    end
    
    if length(Yn_u{g,2}) ~= N
    fs = 1;
    [Yn_u{g,4}, Yn_u{g,5}] = resample(Yn_u{g,1},Yn_u{g,3},fs);
    
    elseif length(Yn_u{g,2}) == N
        Yn_u{g,4} = Yn_u{g,1};
        Yn_u{g,5} = Yn_u{g,3};
    end
    
        if length(Yn_u{g,5}) ~= N
            n = Yn_u{g,5}(1);
            add_begin = linspace(0,n-1,n);   

            if isempty(add_begin) == 0 
                Yn_u{g,5} = [add_begin                   Yn_u{g,5}];
                Yn_u{g,4} = [zeros(1,length(add_begin))  Yn_u{g,4}];
            end

            % m = Yn_u{1,3}(end) - Yn_u{g,5}(end);
            m = (length(time)-1) - Yn_u{g,5}(end);
            if m < 0 || m == 0
                Yn_u{g,5} = Yn_u{g,5}(1:(end-abs(m)));
                Yn_u{g,4} = Yn_u{g,4}(1:(end-abs(m)));
            else
                add_end = linspace(1+Yn_u{g,5}(end), N - 1, m);
                if isempty(add_end) == 0
                    Yn_u{g,5} = [Yn_u{g,5}  add_end];
                    Yn_u{g,4} = [Yn_u{g,4}  zeros(1,length(add_end))];
                end
            end
        end
    
    fclose(fileID);
    g = g + 1; % loop counter
    
catch
   Yn_u{g,4} = zeros(1,N);
   Yn_u{g,4} = zeros(1,N);
   
   fclose(fileID);
   g = g + 1; % loop counter
end
    
end
cd('C:\Users\Asus\Desktop\WRF_BA_5'); % Go back to home dir