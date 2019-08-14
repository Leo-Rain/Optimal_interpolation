function a = ridge(Rjj, Rij,   Xj_u,           Xj_v,           Yn_u,       Yn_v,        beta_u, beta_v)
%        a = ridge(Rjj, Rij_i, Xj_u_Month_wrf, Xj_v_Month_wrf, Yn_u_Sep_c, Yn_v_Sep_c,  beta);

Diff_u = Yn_u - Xj_u;
Diff_v = Yn_v - Xj_v;

Var_Diff_u = var(Diff_u,0,2); % var along the second dim with default normalization 
Var_Diff_v = var(Diff_v,0,2); % var along the second dim with default normalization 

Sigma_Var_Diff_u = diag(Var_Diff_u);
Sigma_Var_Diff_v = diag(Var_Diff_v);

% Sigma_average = mean([Sigma_average_u Sigma_average_v]);
% Sigma_average = sqrt(Sigma_average);

% beta = 5; % CHOICE OF A HYPER PARAMETER!

Sigma_u = Sigma_Var_Diff_u.*beta_u;
Sigma_v = Sigma_Var_Diff_v.*beta_v;

size_Rjj = size (Rjj);

Zero_matrix = zeros(size_Rjj(1)/2,size_Rjj(1)/2);

Aug_matrix = [Sigma_u Zero_matrix; Zero_matrix Sigma_v];

% Aug_matrix = eye(size_Rjj).*beta;

Rjj_ridge_add = Rjj + Aug_matrix;

Size_Rij = size (Rij);

a = zeros(Size_Rij(1),Size_Rij(2),Size_Rij(3));
for y = 1:Size_Rij(1)
    for x = 1:Size_Rij(2)
        Rij_yx = squeeze(Rij(y,x,:))';
        a_yx = Rij_yx/Rjj_ridge_add; % Resolving the linear system !!!!

        a(y,x,:) = a_yx;
%       a(y,x,:) = Rij_yx * pinv(Rjj_ridge_add);
        
    end
end

end