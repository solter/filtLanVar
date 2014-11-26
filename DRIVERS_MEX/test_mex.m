%% a test driver for MEX code "laneig" and "filtlan"


% load the (sparse) matrix
disp('load ../DATA/laplacian221.mat');
load ../DATA/laplacian221.mat


%% the Lanczos algorithm (w/ partial reorthogonalization) for eigenvalue computations

% compute the 7 smallest eigenvalues
% all less than 0.5
fprintf('\n[V1, D1, info1] = laneig(A, 7, ''sa'');\n');
[V1, D1, info1] = laneig(A, 7, 'sa');
D1
fprintf('||A*V1-V1*D1||_F = %e\n\n', norm(A*V1-V1*D1,'fro'));

% computes the 17 smallest eigenvalues
% the first 7 less than 0.5, the next 9 less than 1.0, the last one greater than 1.0
fprintf('\n[V2, D2, info2] = laneig(A, 17, ''sa'');\n');
[V2, D2, info2] = laneig(A, 17, 'sa');
D2
fprintf('||A*V2-V2*D2||_F = %e\n\n', norm(A*V2-V2*D2,'fro'));


%% the filtered Lanczos algorithm (w/ partial reorthogonalization) for eigenvalue computations

% compute the eigenvalues in (-inf,0.5], with polydeg==10
% in total there are 7 eigenvalues
fprintf('\n[V3, D3, info3] = filtlan(A, [-inf,0.5], 10);\n');
[V3, D3, info3] = filtlan(A, [-inf,0.5], 10);
D3
fprintf('||A*V3-V3*D3||_F = %e\n\n', norm(A*V3-V3*D3,'fro'));

% compute the eigenvalues in (-inf,1.0], with polydeg==10
% in total there are 16 eigenvalues
fprintf('\n[V4, D4, info4] = filtlan(A, [-inf,1.0], 10);\n');
[V4, D4, info4] = filtlan(A, [-inf,1.0], 10);
D4
fprintf('||A*V4-V4*D4||_F = %e\n\n', norm(A*V4-V4*D4,'fro'));

% compute the eigenvalues in [0.5,1.0], with pollydeg==20
% in total there are 9 eigenvalues
fprintf('\n[V5, D5, info5] = filtlan(A, [0.5,1.0], 20);\n');
[V5, D5, info5] = filtlan(A, [0.5,1.0], 20);
D5
fprintf('||A*V5-V5*D5||_F = %e\n\n', norm(A*V5-V5*D5,'fro'));
