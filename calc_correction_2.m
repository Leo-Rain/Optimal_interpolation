
function [d_u, d_v] = calc_correction_2(a, Xj_u_wrf, Xj_v_wrf, Yn_u_c, Yn_v_c, X, Y, M, T, N)

%%%%% Spare vector %%%%%%%%%%%%%
% N = length(Xj_u_wrf);
% step = 5;
% Xj_u_wrf = Xj_u_wrf(:,1:step:N);
% Xj_v_wrf = Xj_v_wrf(:,1:step:N);
% 
% Yn_u_c = Yn_u_c(:,1:step:N);
% Yn_v_c = Yn_v_c(:,1:step:N);
% 
% N = length(Xj_u_wrf); % changed length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
Diff_shift_u = zeros(M, N, (2*T+1));
Diff_shift_v = zeros(M, N, (2*T+1));

for Shift = -T:T
    
%     Diff_shift_u(:,:,Shift+T+1) =  Yn_u_c - circshift(Xj_u_wrf,Shift,2);
%     Diff_shift_v(:,:,Shift+T+1) =  Yn_v_c - circshift(Xj_v_wrf,Shift,2);
%     
Diff_shift_u(:,:,Shift+T+1) =  circshift(Yn_u_c,Shift,2) - circshift(Xj_u_wrf,Shift,2);
Diff_shift_v(:,:,Shift+T+1) =  circshift(Yn_v_c,Shift,2) - circshift(Xj_v_wrf,Shift,2);
    
end
    
%   A1 = squeeze(Diff_shift_u(1,:,:));
%   A2 = squeeze(Diff_shift_u(2,:,:));
%   
%   A3 = squeeze(Diff_shift_v(1,:,:));
%   A4 = squeeze(Diff_shift_v(2,:,:));
    
Size_Diff = size(Diff_shift_u);

if T ~= 0
    S3 = Size_Diff(3);
else
    S3 = 1;
end

S1 = Size_Diff(1);
    
Diff_shift_u_cat = zeros( N,(S1*S3) );
Diff_shift_v_cat = zeros( N,(S1*S3) );

for i = 1:S1
    Diff_shift_u_cat(:, (1+S3*(i-1)): (S3+S3*(i-1)) ) = squeeze(Diff_shift_u(i,:,:));
    Diff_shift_v_cat(:, (1+S3*(i-1)): (S3+S3*(i-1)) ) = squeeze(Diff_shift_v(i,:,:));
end

% Split a in two: a_u and a_v:
Size_a = size(a);  % Size_a = [439 479 84]
la = Size_a(3);    % la = 84
half = ceil(la/2); % for odd number of bit-stream length (half = 42)

a_u = a(:,:,1:half);
a_v = a(:,:,half + 1 : end);

d_u = zeros(Y,X,N);
d_v = zeros(Y,X,N);

for y = 1:Y
    for x = 1:X
        a_u_yx = squeeze(a_u(y,x,:));
        d_u(y,x,:) = Diff_shift_u_cat * a_u_yx;
        
        a_v_yx = squeeze(a_v(y,x,:));
        d_v(y,x,:) = Diff_shift_v_cat * a_v_yx; 
        
        if y == 300 && x == 300
            foo = 0;
        end
        
    end
end


