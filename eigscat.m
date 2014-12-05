% Matlab code eigscat.m
% For "Applied Numerical Linear Algebra",  Question 4.14
% Written by James Demmel, Oct 25, 1995
%                Modified, Jun  5, 1997
%                Modified, Oct 15, 2004
%
%  Plot eigenvalues of a matrix and eigenvalues of
%  random perturbations of the same matrix. 
%  This illustrates the sensitivity of eigenvalues to perturbations.
%
%  Inputs:
%     a = input matrix
%     err = size of perturbation
%     m = number of perturbed matrices to compute
%
%  Outputs:
%     Three plots where each symbol is the location of an eigenvalue of 
%     a perturbed matrix:
%        'o' marks the location of each unperturbed eigenvalue. 
%        'x' marks the location of each perturbed eigenvalue, where a
%                  real perturbation matrix of norm err is added to a. 
%        '.' marks the location of each perturbed eigenvalue, where a
%                  complex perturbation  matrix of norm err is added to a. 
%

%  Compute data
%  ea contains eigenvalues of unperturbed matrix a
%  er contains eigenvalues of matrix with real perturbations
%  ei contains eigenvalues of matrix with complex perturbations

function eigscat(a, err, m, radius)
    ea=conj(eig(a)'); % Why the conj and the transpose? 
    er=ea; 
    ei=ea;
    
    for i=1:m
        r=rand(size(a))-.5; 
        r=(err/norm(r))*r; 
        er=[er,conj(eig(a+r))']; 
    end
    j=sqrt(-1); 
    for i=1:m
        r=(rand(size(a))-.5)+j*(rand(size(a))-.5); 
        r=(err/norm(r))*r;
        ei=[ei,conj(eig(a+r)')]; 
    end
    %  Plot data
    axis('square'),
    % Real perturbations
    hold off,plot(real(er),imag(er),'xr'),hold on,grid,
    plot(real(ea),imag(ea),'ok'),title('Real Perturbations'),
    for i = 1:length(a)
        hold on, circle(real(ea(i)),imag(ea(i)), radius(i)), shg;
    end
    disp('pause (hit return to continue)'),pause
    % Complex perturbations
    hold off,plot(real(ei),imag(ei),'.b'),hold on,grid,
    plot(real(ea),imag(ea),'ok'),title('Complex Perturbations'),
    for i = 1:length(a)
        hold on, circle(real(ea(i)),imag(ea(i)), radius(i)), shg;
    end
    disp('pause (hit return to continue)'),pause
    % Both types of perturbation
    hold off,plot(real(ei),imag(ei),'.b'),hold on,grid,
    plot(real(er),imag(er),'xr'),plot(real(ea),imag(ea),'ok'),
    title('Real and Complex Perturbations')
    for i = 1:length(a)
        hold on, circle(real(ea(i)),imag(ea(i)), radius(i)), shg;
    end
    %
    %  Compute and print condition numbers
    [n,n]=size(a); cnd=1; format short e,
    [v,d]=eig(a);for i=1:n, v(:,i)=v(:,i)/norm(v(:,i)); end
    vi=inv(v); for i=1:n, cnd(i)=norm(vi(i,:)); end
    disp('       Eigenvalues      Condition Numbers'), [diag(d),cnd']
end