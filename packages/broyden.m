function [xv,it]=broyden(x,f,n,tol)
% Broyden's method for solving a system of n non-linear equations
% in n variables.
%
% Example call: [xv,it]=broyden(x,f,n,tol)
% Requires an initial approximation column vector x. tol is required
% accuracy. User must define function f, for example see page 115.
% xv is the solution vector, parameter it is number of iterations
% taken. WARNING. Method may fail, for example, if initial estimates
% are poor.
%
x=x';
fr=zeros(n,1); it=0; xv=x;
%Set initial Br
Br=eye(n);
fr=feval(f, xv);
fr=fr';
while norm(fr)>tol
    it=it+1;
    if it>100
        it
        error('it large');
    end
    pr=-Br*fr;
    tau=1;
    xv1=xv+tau*pr;
    xv=xv1;
    oldfr=fr;
    fr=feval(f,xv);
    fr=fr';
    %Update approximation to Jacobian using Broydens formula
    y=fr-oldfr;
    oldBr=Br;
    oyp=oldBr*y-pr; pB=pr'*oldBr;
    M=oyp*pB;
    Br=oldBr-M./(pr'*oldBr*y);
end;
it

