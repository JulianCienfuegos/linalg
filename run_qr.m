dim = 30;
A = rand(dim); A = A'*A;
m=1000;
[A,err] = qrplt(A,m);
plot(1:(m-1), log(err))
length(err)