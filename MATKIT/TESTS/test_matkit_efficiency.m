% This m-file performs similar operations to those in
% "test_matkit_efficiency.cpp", which is for testing the efficiency of some
% MATKIT functions.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    efficiency test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('### efficiency test\n\n');

% sparse matrix - vector multiplication
nr = 20000;
nc = 30000;
bnd = 50;
S = sparse(nr, nc);
for j = 1:nc
    for i = max(1,j-bnd):min(j+bnd,nr)
        S(i,j) = 1/(1.0+abs(i-j));
    end
end
fprintf('S is of size %d-by-%d with nnz = %d\n', size(S,1), size(S,2), nnz(S));
fprintf('||S||_1   = %g\n', norm(S,1));
fprintf('||S||_F   = %g\n', norm(S,'fro'));
fprintf('||S||_inf = %g\n\n', norm(S,inf));

% sparse matrix - vector multiplication
fprintf('v = ones(%d, 1);\n', nr);
v = ones(nc, 1);
fprintf('w=S*v time used (1000 times): ');
t = cputime;
for i = 1:1000
    w = S*v;
end
fprintf('%g seconds\n', cputime-t);
fprintf('||w||_1   = %g\n', norm(w,1));
fprintf('||w||_2   = %g\n', norm(w,2));
fprintf('||w||_inf = %g\n\n', norm(w,inf));

% vector - sparse matrix multiplication
fprintf('v2 = ones(1, %d);\n', nc);
v2 = ones(1, nr);
fprintf('w2=v2*S time used (1000 times): ');
t = cputime;
for i = 1:1000
    w2 = v2*S;
end
fprintf('%g seconds\n', cputime-t);
fprintf('||w2||_1   = %g\n', norm(w2,1));
fprintf('||w2||_2   = %g\n', norm(w2,2));
fprintf('||w2||_inf = %g\n\n', norm(w2,inf));

% dense matrix transpose
dim = 10000;
M = rand(dim, dim);
fprintf('time used to transpose a %d-by-%d matrix: ', dim, dim);
t = cputime;
transpose(M);
fprintf('%g seconds\n', cputime-t);

% dense matrix - vector multiplication
v = rand(dim, 1);
fprintf('dense matrix - vector multiplication (n=%d) time used: ', dim);
t = cputime;
w = M*v;
fprintf('%g seconds\n', cputime-t);

v = rand(1, dim);
fprintf('vector - dense matrix multiplication (n=%d) time used: ', dim);
t = cputime;
w = v*M;
fprintf('%g seconds\n', cputime-t);

% dense matrix - dense matrix multiplication
dim = 2000;
M = rand(dim);
fprintf('dense matrix - dense matrix multiplication (n=%d) time used: ', dim);
t = cputime;
N = M*M;
fprintf('%g seconds\n', cputime-t);
