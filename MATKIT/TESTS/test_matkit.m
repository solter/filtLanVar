% This m-file performs similar operations to those in "test_matkit.cpp",
% which illustrates how to use MATKIT classes Vector, Matrix,
% SymmetricMatrix, and SparseMatrix.
% For experienced OCTAVE/MATLAB users, comparing this m-file with
% test_matkit.cpp may help understand how to use MATKIT.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    machine precision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('### machine precision\n\n');

fprintf('infinity  = %g\n', Inf);
fprintf('-infinity = %g\n', -Inf);
fprintf('eps = %g\n', eps);
fprintf('1 + eps   - 1 = %g\n', 1+eps-1);
fprintf('1 + eps/2 - 1 = %g\n', 1+eps/2-1);

fprintf('\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    vector routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('### vector routines\n\n');

% vector v
v = [ 1;  2;  3 ];
v

% subvector of v
fprintf('v(3:-1:1)\n');
v(3:-1:1)
fprintf('v(2:3)\n');
v(2:3)

% ||v||^2
fprintf('||v||^2 = %g\n\n', norm(v)^2);

% diag(v), etc.
fprintf('diag(v)\n');
diag(v)
fprintf('diag(v,1)\n');
diag(v,1)
fprintf('sparse(diag(v,-2))\n')
sparse(diag(v,-2))

% vector - scalar operations
fprintf('-v\n');
-v
fprintf('3*v\n');
3*v
fprintf('v/2\n');
v/2
fprintf('1+v\n');
1+v
fprintf('v-2\n');
v-2
fprintf('v = v*3;\n');
v = v*3;
v
fprintf('v = v-2;\n');
v = v-2;
v
fprintf('v = v+5;\n');
v=v+5;
v
fprintf('v = v/3;\n');
v=v/3;
v

% vector - vector addition/subtraction
w = [ 1; -2; 4 ];
w
fprintf('v+w\n');
v+w
fprintf('v-w\n');
v-w
fprintf('v = v-w;\n');
v=v-w;
v
fprintf('v = v+w;\n');
v=v+w;
v

% vector inner/outer products
fprintf('v''*w\n');
v'*w
fprintf('v*v''\n');
v*v'
fprintf('v*w''\n');
v*w'
fprintf('v*w''+w*v''\n');
v*w'+w*v'

fprintf('\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    matrix routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('### matrix routines\n\n');

% eye
fprintf('eye(2)\n');
eye(2)
fprintf('eye(3,2)\n');
eye(3,2)

% zeros
fprintf('zeros(2)\n');
zeros(2)
fprintf('zeros(3,2)\n');
zeros(3,2)

% ones
fprintf('ones(3)\n');
ones(3)
fprintf('ones(3,1)\n');
ones(3,1)

% sparse(ones(2,3)-eye(2,3))
fprintf('sparse(ones(2,3)-eye(2,3))\n');
sparse(ones(2,3)-eye(2,3))

% A34(i,j) = (i-1)*4 + j for i = 1,...,3 and j = 1,...,4
A34 = [ 1  2  3  4;  5  6  7  8;  9  10  11  12 ];
A34

% submatrix of A34
fprintf('A34(:,4)\n');
A34(:,4)
fprintf('A34(2,:)\n');
A34(2,:)
fprintf('A34(:,4:-1:2)\n');
A34(:,4:-1:2)
fprintf('A34(:,1:3)\n');
A34(:,1:3)
fprintf('A34(:,[3,1])\n');
A34(:,[3,1])
fprintf('A34(2:3,:)\n');
A34(2:3,:)
fprintf('A34(3:-1:1,:)\n');
A34(3:-1:1,:)
fprintf('A34([3,1],:)\n');
A34([3,1],:)
fprintf('A34(3:-1:2,3:-1:1)\n');
A34(3:-1:2,3:-1:1)
fprintf('A34(1:2,2:4)\n');
A34(1:2,2:4)

% diag of A34
fprintf('diag(A34)\n');
diag(A34)
fprintf('diag(A34,2)\n');
diag(A34,2)
fprintf('diag(A34,-1)\n');
diag(A34,-1)

% transpose of A34
fprintf('transpose(A34)\n');
transpose(A34)

% norms of A34
fprintf('||A34||_F^2 = %g\n', norm(A34,'fro')^2);
fprintf('||A34||_1   = %g\n', norm(A34,1));
fprintf('||A34||_inf = %g\n\n', norm(A34,inf));

% miscellaneous operations
fprintf('A24 = A34(1:2,:);\n')
A24 = A34(1:2,:);
A24
fprintf('A42 = transpose(A24);\n');
A42 = A24';
A42
fprintf('[A42,zeros(4,2)] + [A24;zeros(2,4)]\n');
[A42,zeros(4,2)] + [A24;zeros(2,4)]

% permutation of rows/columns of A34
perm3 = [ 2  3  1 ];
perm3
perm4 = [ 3  2  4  1 ];
perm4
fprintf('Z(:,perm4) = A34;\n');
Z(:,perm4) = A34;
Z
fprintf('Z(perm3,:) = A34;\n');
Z(perm3,:) = A34;
Z
fprintf('Z(perm3,perm4) = A34;\n');
Z(perm3,perm4) = A34;
Z

% matrix - matrix multiplication
fprintf('A22 = A24*A42;\n');
A22 = A24*A42;
A22
fprintf('A22 = A22-50;\n');
A22 = A22-50;
A22
fprintf('A22 = A22*0.5;\n');
A22 = A22*0.5;
A22
fprintf('A22 = A22/2;\n');
A22 = A22/2.0;
A22
fprintf('A22 = A22+1;\n');
A22 = A22+1.0;
A22
fprintf('A24b = A22*A24;\n');
A24b = A22*A24;
A24b
fprintf('A24b - A24\n');
A24b-A24
fprintf('A24b = A24b+A24;\n');
A24b = A24b+A24;
A24b
fprintf('[A24,A22]\n');
[A24,A22]
fprintf('[A42;A22]\n');
[A42;A22]

fprintf('\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    matrix - vector operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('### matrix - vector operations\n\n');

% A(i,j) = i+(j-1)*3 for i,j = 1,...,3
A = [ 1  4  7;  2  5  8;  3  6  9 ];

% matrix - vector multiplication
A
v
fprintf('A*v\n');
A*v
fprintf('v''*A\n');
v'*A

fprintf('\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    symmetric matrix routines (simulated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('### symmetric matrix routines (simulated)\n\n');

B = [ 1  2  4;  2  3  5;  4  5  6 ];
B

% submatrix of B
fprintf('B(:,3)\n');
B(:,3)
fprintf('B(2,:)\n');
B(2,:)
fprintf('B(:,2:3)\n');
B(:,2:3)
fprintf('B(2:-1:1,:)\n');
B(2:-1:1,:)
fprintf('B(2:3,2:-1:1)\n');
B(2:3,2:-1:1)
fprintf('B(2:3,2:3)\n');
B(2:3,2:3)
fprintf('B(2:-1:1,2:-1:1)\n');
B(2:-1:1,2:-1:1)

% diag of B
fprintf('diag(B)\n');
diag(B)
fprintf('diag(B,1)\n');
diag(B,1)
fprintf('diag(B,-2)\n');
diag(B,-2)
fprintf('\n');

% transpose of B
fprintf('transpose(B)\n');
transpose(B)

% norms of B
fprintf('||B||_F^2 = %g\n', norm(B,'fro')^2);
fprintf('||B||_1   = %g\n', norm(B,1));
fprintf('||B||_inf = %g\n\n', norm(B,inf));

% symmetric permutation of B
perm3
fprintf('C(perm3,perm3) = B;\n');
C(perm3,perm3) = B;
C
fprintf('B+C\n');
B+C
fprintf('B-C\n');
B-C
fprintf('B*C\n');
B*C

fprintf('\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    symmetric matrix - vector operations (simulated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('### symmetric matrix - vector operations (simulated)\n\n');

v
fprintf('B*v\n');
B*v
fprintf('B''*v\n');
v'*B

fprintf('\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    matrix - symmetric matrix operations (simulated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('### matrix - symmetric matrix operations (simulated)\n\n');

% matrix - symmetric matrix addition/subtraction
A
fprintf('A+B\n');
A+B
fprintf('A = A-B;\n');
A = A-B;
A
fprintf('B+A\n');
B+A
fprintf('B-A\n');
B-A

% preparation for matrix - symmetric matrix multiplication
A23 = [ 1  4  7;  2  6  8 ];
A23
fprintf('');
A32 = transpose(A23);
A32

% matrix - symmetric multiplication
fprintf('B*A32\n');
B*A32
fprintf('A23 = A23*B;\n');
A23 = A23*B;
A23

fprintf('\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    sparse matrix routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('### sparse matrix routines\n\n');

% sparse matrix creation and element access
S22 = 3*speye(2);
S22
fprintf('S22(1,2) == %g\n', full(S22(1,2)));
fprintf('S22(2,2) == %g\n\n', full(S22(2,2)));

S32 = speye(3,2)/2;
S32
fprintf('S32(2,2) = S32(2,2)*2;\n');
S32(2,2) = S32(2,2)*2;
fprintf('S32(3,1) = 5.0;\n');
S32(3,1) = 5;
S32

% convert S32 to a dense matrix
fprintf('full(S32)\n');
full(S32)

% Frobenius norm of S32
fprintf('||S32||_F^2 = %g\n\n', norm(S32,'fro')^2);

% diag of S32
fprintf('full(diag(S32))\n');
full(diag(S32))
fprintf('full(diag(S32,1))\n');
full(diag(S32,1))
fprintf('full(diag(S32,-2))\n');
full(diag(S32,-2))
fprintf('\n');

% transpose of S32, etc.
fprintf('S32(3,2) = -1;\n');
S32(3,2) = -1;
S32
fprintf('transpose(S32)\n');
transpose(S32)

% permutation of rows and/or columns of S32
perm2 = [ 2  1 ];
perm2
perm3
fprintf('R32(:,perm2) = S32;\n');
R32(:,perm2) = S32;
R32
fprintf('R32(perm3,:) = S32;\n');
R32(perm3,:) = S32;
R32
fprintf('R32(perm3,perm2) = S32;\n');
R32(perm3,perm2) = S32;
R32

% sparse matrix - sparse matrix addition/subtraction
T32 = sparse([ 0  0.5;  0.2  2.1;  0  0 ]);
T32
fprintf('U32 = S32 + T32;\n');
U32 = S32 + T32;
U32
fprintf('T32 = T32 - S32;\n');
T32 = T32 - S32;
T32

% submatrix of U32
fprintf('U32(:,2:-1:1)\n');
U32(:,2:-1:1)
fprintf('U32(3:-1:2,:)\n');
U32(3:-1:2,:)
fprintf('U32(3:-1:2,2:-1:1)\n');
U32(3:-1:2,2:-1:1)

% construction of sparse tridiagonal matrix
fprintf('T = sparse(diag([2 3],-1) + diag([3 4 5]) + diag([1 2],1));\n');
T = sparse(diag([2 3],-1) + diag([3 4 5]) + diag([1 2],1));
T
fprintf('\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    sparse matrix - vector operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('### sparse matrix - vector operations\n\n');

v = [ 1;  2 ];
v
fprintf('w = S32*v;');
w = S32*v;
w
fprintf('w = w''*S32');
w = w'*S32;
w

fprintf('\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    matrix - sparse matrix operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('### matrix - sparse matrix operations\n\n');

% matrix - sparse matrix multiplication
S32
A23
fprintf('S32*A23\n');
S32*A23
fprintf('A23*S32\n');
A23*S32

fprintf('\n\n');
