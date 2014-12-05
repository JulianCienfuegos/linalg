dim = 50;
A = rand(dim); A = A'*A;
m=10000;
[A,err] = qrplt(A,m);
err
plot(1:m, log(err))