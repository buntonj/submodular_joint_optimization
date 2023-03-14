%% FAST NON-NEGATIVE QP SOLVING
% Solves optimization problems of the form:
%
%  minimize  0.5*x'*Q*x + p*x
%    s.t.       x >= 0
%
% inputs:
% Q - n x n positive semidefinite matrix
% p - n x 1 array
% tol - tolerance on KKT conditions for stopping criteria
%
% returns:
% x - n x 1 array solution of optimization problem
% val - scalar, achieved cost of solution

function [x,val] = fnnqp(Q,p,tol)
if size(Q) == 0
    val = 0;
    x = [];
    return
end
n = size(p,1);
x = zeros(n,1);
I = false(n,1); %inactive constraints (i.e. x > 0)
A = true(n,1); %active constraints  (i.e. x=0)
w  = -p - Q*x;
s = zeros(n,1);
[m,j] = max(w);
iter = 0;
iter_max = 30*n;
while (m > tol) & any(A,'all')
    I(j) = true;
    A(j) = false;
    s(I) = -Q(I,I)\ p(I);
    s(A) = 0;
    while min(s(I)) <= 0 && iter < iter_max
        iter = iter + 1;
        R = (s <= tol) & I;
        alpha = min(x(R)./(x(R)-s(R)));
        x = x + alpha*(s-x);
        I = x>0;
        A = ~I;
        s(I) = -Q(I,I)\ p(I);
        s(A) = 0;
    end
    x = s;
    w = -p - Q*x;
    [m,j] = max(w);
end
val = 0.5*x'*Q*x + p'*x;
end