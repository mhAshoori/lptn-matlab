function [x] = SOR_HW(A,b,x_0,omega)
% Input a square matrix A, b, initial x and value of omega
format long;
N = 1000; %number of iteration
n = length(A);
tol = 0.0001;
x =zeros(n,1);
%Decomposing the Square matrix A into three matrices: diagonal matrix (D); strictly lower triangular matrix (L); strictly upper triangular matrix(U)
D = diag(diag(A));
L =-tril(A,-1);
U = -triu(A,1);
a = (D-omega*L);
for i=1:N
    x = a\(((1-omega)*D + omega*U)*x_0) + omega*(a\b);
    if norm(x-x_0)<tol
        break;
    end
    x_0=x;
end
end