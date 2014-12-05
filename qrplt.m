% Matlab code qrplt.m
% For "Applied Numerical Linear Algebra",  Question 4.15
% Written by James Demmel, Oct 25, 1995
%                Modified, Jun  5, 1997
%
%  Plot the diagonal entries of a matrix undergoing unshifted QR iteration
%
%  Inputs:
%    a = matrix
%    m = number of QR iterations
%
%  Outputs:
%    n curves, one for each diagonal entry of a 
%
function [a, err] = qrplt(a, m)
    hold off
    e=diag(a);
    for i=1:m,
       [q,r]=qr(a);dd=diag(sign(diag(r)));r=dd*r;q=q*dd;a=r*q; ...
       e=[e,diag(a)];
    end
    clf
    plot(e','k'),grid
    title('plot of each diagonal matrix entry during QR iteration')
    shg
    [~,m] = size(e);
    err = zeros((m-1), 1); 
    for col = 1:(m-1)
        err(col) = norm(e(:,m) - e(:,col))/norm(e(:,m));
    end
end