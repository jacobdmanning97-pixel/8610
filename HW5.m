%Problem 5a 
rng('default');
 F = randn(10,10);

K=HBD(F);
Sig=iterGKbdiag(K);
[U,S,V] = svd(K);

%Problem 5b
col = linspace(-1,1,1024*1024+1)';
A = col.^(0:31);

R=HQR(A);
B=HBD(R);

[m,n]=size(B);

block = [zeros(m,n) B';B zeros(n,m)];
eblock = eig(block);
eblock = sort(eblock);
eblock = eblock(n+1:2*n);
sblock = sqrt([eblock(1:5);eblock(n-4:n)]);

e = eig(A'*A);
e = sort(e);
s = sqrt([e(1:5);e(n-4:n)]);

function R = HQR(A)
    [m,n] = size(A);
    R = A;

    for k = 1:n
        x = R(k:m,k);
        e = zeros(length(x),1);
        e(1) = norm(x);
        if x(1) == 0
            beta = 1;
        else
            beta = sign(x(1));
        end
        u = beta*e + x;
        v = u / norm(u);
        R(k:m,k:n) = R(k:m,k:n) - 2*v*(v'*R(k:m,k:n));
    end

    R = R(1:n,1:n);
end

function R = HBD(A)
    [m,n] = size(A);
    R = A;

    for k = 1:n
        %Left Hand Side
        x = R(k:m,k);
        e = zeros(length(x),1);
        e(1) = norm(x);
        if x(1) == 0
            beta = 1;
        else
            beta = sign(x(1));
        end
        u = beta*e + x;
        v = u / norm(u);
        R(k:m,k:n) = R(k:m,k:n) - 2*v*(v'*R(k:m,k:n));

        %Right Hand Side
        if k<n-1
            x = R(k,k+1:n);
            e = zeros(1,length(x));
            e(1) = norm(x);
            if x(1) == 0
                beta = 1;
            else
                beta = sign(x(1));
            end
            u = beta*e + x;
            v = u / norm(u);
            R(k:m,k+1:n) = R(k:m,k+1:n) - 2*R(k:m,k+1:n)*v'*v;
        end
    end

    R = R(1:n,1:n);
       
    %Clean to use GKbdiag code (does not recognize e_mach as zero)
    for i = 1:m
        for j = 1:n
            if abs(R(i,j))<10^-15
                R(i,j)=0;
            end
        end
    end
end

function svals = iterGKbdiag(B)
% an illustrative code for computing the singular values of a bidiagonal B
nsvd = 0;
n = min(size(B));
bsize = n;
svals = zeros(n,1);
    while nsvd < n
        if bsize >= 2
            B = onestepGKbdiag(B);
            if abs(B(end-1,end)) < eps/2*(abs(B(end-1,end-1))+abs(B(end,end)))
                nsvd = nsvd + 1;
                svals(nsvd) = abs(B(end,end));
                B = B(1:end-1,1:end-1);
                bsize = bsize - 1;
            end
        else
            nsvd = nsvd + 1;
            svals(nsvd) = abs(B(1,1));
        end
    end
end

function [D,CSL,CSR] = onestepGKbdiag(B)
% following Golub & van Loan Alg. 8.6.1: 
% apply Givens rotation G' on the left (pre-multiply by [c -s; s c]), and 
% apply G on the right (post-multiply by [c s; -s c]). 
n = min(size(B));
if nnz(tril(B,-1)) > 0 || nnz(triu(B,2)) > 0
    error('Input matrix B is not upper bidiagonal');
end
if n >= 3
    T = [B(n-1,n-1)^2+B(n-2,n-1)^2 B(n-1,n-1)*B(n-1,n); 
        B(n-1,n-1)*B(n-1,n) B(n,n)^2+B(n-1,n)^2];
else
    T = B'*B;
end
CSL = zeros(2,n-1); CSR = zeros(2,n-1);
evT = eig(T);
[~,idx] = min(abs(evT-(B(n,n)^2+B(n-1,n)^2)));
lambda = evT(idx);
u = [B(1,1)^2-lambda; B(1,1)*B(1,2)];
[c,s] = givens(u);
B(1:2,1:2) = B(1:2,1:2)*[c s; -s c];
CSR(:,1) = [c; s];

for k = 1 : n-2
    [c,s] = givens(B(k:k+1,k));
    B(k:k+1,k:k+2) = [c s; -s c]'*B(k:k+1,k:k+2);
    B(k+1,k) = 0;
    CSL(:,k) = [c; s];

    [c,s] = givens(B(k,k+1:k+2)');
    B(k:k+2,k+1:k+2) = B(k:k+2,k+1:k+2)*[c s; -s c];
    B(k,k+2) = 0;
    CSR(:,k+1) = [c; s];
end
[c,s] = givens(B(n-1:n,n-1));
B(n-1:n,n-1:n) = [c -s; s c]*B(n-1:n,n-1:n);
B(n,n-1) = 0;
D = B;

end

function [c,s] = givens(u)
    if u(2) == 0
        c = 1;  s = 0;
    else
        if abs(u(2)) > abs(u(1))
            tau = -u(1)/u(2);   s = 1/sqrt(1+tau^2);    c = s*tau;
        else
            tau = -u(2)/u(1);   c = 1/sqrt(1+tau^2);    s = c*tau;
        end
    end
end