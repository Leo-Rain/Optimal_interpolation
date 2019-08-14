function Xj = corrsp_field_serie(A,M,N,Xj_c)
%-------------------------------------------------------------------------
% This function extracts from WRF the time-series corresponding to RP5
% stations coordinates
%-------------------------------------------------------------------------
Xj = zeros(M,N);      % M points, each N values long
%Xj_cent = zeros(M,N); % M points, each N values long
for i = 1:M 
          Xj(i,:) = reshape( A(Xj_c(i,1), Xj_c(i,2),:) ,1,[]); % Make a linear array (?): WHEN EXTRACT FROM A - SWAP X and Y!!!!!!
   % Xj_cent(i,:) = reshape( A_cent(Xj_c(i,1),     Xj_c(i,2),:) ,1,[]); % Make a linear array (for xcorr)
end
end