function A = laplace(N)
%
% Diskreter Laplace-Operator (-\Delta), 2 Dimensionen
%
e = ones(N,1);
E = spdiags(e,0,N,N);
F = spdiags([-1*e 2*e -1*e],[-1 0 1],N,N);

A = kron(F,E)+kron(E,F);

A=full(A);