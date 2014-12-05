close all; clear all;
which_example = 1

%{
Theorem 4.4 offers a bound on the perturbation in lambda given a
perturbation of A.
 
I am neglecting the second order term in the bound, because the 
perturbations are so small.
%}
n = 5;
err = 1e-8;
m = 100;
radius = zeros(n,1);
quit = 0;

if(which_example==1)    
    A = 100*rand(n); % A random, nondefective matrix
    [R,D1] = eig(A); R = R/norm(R); % Right eigenvectors
    [L,D2] = eig(A.'); L = conj(L)/norm(conj(L)); % Left eigenvectors
    for i = (1:n)
        radius(i) = err / abs(L(:,i)'*R(:,i));
        if(radius(i) > err*100)
            quit = 1;
        end
    end
    if(quit == 0)
        disp('radii:')
        disp(radius);
        eigscat(A, err,m, radius);
    end
    
end

if(which_example==2)
    A = 100*rand(n);
    A(:,1) = 2*A(:,2); % Make A defective
    [R,D1] = eig(A); R = R/norm(R); % Right eigenvectors
    [L,D2] = eig(A.'); L = conj(L)/norm(conj(L)); % Left eigenvectors
    for i = (1:n)
        radius(i) = err / abs(L(:,i)'*R(:,i));
        if(radius(i) > err*100)
            quit = 1;
        end
    end
    if(quit == 0)
        disp('radii:')
        disp(radius);
        eigscat(A, err,m, radius);
    end
    
end
