load west0479;
C=full(west0479);

[H,Q,iter] = HW8_QReig(C);
E=ordeig(H);
E=sort(E);
D=eig(C);
D=sort(D);
norm(E-D)/norm(D);

u=cos((0:2048)/2048*pi);
B=vander(u);
flag=true;

[V,H]=AR(B,50,flag);

norm(V'*V-eye(51))

l=zeros(11,4);
flag=true;

m=30;
[V,H]=AR(A,m,flag);
[Veig,Deig]=eig(H(1:end-1,:));
[Deig,idx]=sort(diag(Deig),'descend');
Veig=Veig(:,idx);
Veig=V(:,1:end-1)*Veig;

for i=1:11
    l(i,1)=norm(A*Veig(:,i)-Deig(i)*Veig(:,i))/norm(A*Veig(:,i));
end

tiledlayout(3,2)
nexttile
plot(real(ev),imag(ev),"o")
nexttile
plot(real(Deig),imag(Deig),"o")

m=60;
[V,H]=AR(A,m,flag);
[Veig,Deig]=eig(H(1:end-1,:));
[Deig,idx]=sort(diag(Deig),'descend');
Veig=Veig(:,idx);
Veig=V(:,1:end-1)*Veig;

for i=1:11
    l(i,2)=norm(A*Veig(:,i)-Deig(i)*Veig(:,i))/norm(A*Veig(:,i));
end

nexttile
plot(real(Deig),imag(Deig),"o")

m=100;
[V,H]=AR(A,m,flag);
[Veig,Deig]=eig(H(1:end-1,:));
[Deig,idx]=sort(diag(Deig),'descend');
Veig=Veig(:,idx);
Veig=V(:,1:end-1)*Veig;

for i=1:11
    l(i,3)=norm(A*Veig(:,i)-Deig(i)*Veig(:,i))/norm(A*Veig(:,i));
end

nexttile
plot(real(Deig),imag(Deig),"o")

m=150;
[V,H]=AR(A,m,flag);
[Veig,Deig]=eig(H(1:end-1,:));
[Deig,idx]=sort(diag(Deig),'descend');
Veig=Veig(:,idx);
Veig=V(:,1:end-1)*Veig;

for i=1:11
    l(i,4)=norm(A*Veig(:,i)-Deig(i)*Veig(:,i))/norm(A*Veig(:,i));
end

nexttile
plot(real(Deig),imag(Deig),"o")
l

function [V,H]=AR(A,m,flag)
    [n,~]=size(A);
    r0=zeros(n,1);
    r0(1)=1;
    H=zeros(m+1,m);
    V=zeros(n,m+1);
    V(:,1)=r0/norm(r0);
    for k=1:m
        w=A*V(:,k);
        for i=1:k
            H(i,k)=V(:,i)'*w;
            w=w-H(i,k)*V(:,i);
        end
        if flag==true
            for i=1:k
                dH(i,k)=V(:,i)'*w;
                w=w-dH(i,k)*V(:,i);
                H(i,k)=H(i,k)+dH(i,k);
            end
        end
        H(k+1,k)=norm(w);
        V(:,k+1)=w/norm(w);
    end
end

function [H,Q,iter] = HW8_QReig(A)

[m,n] = size(A);
if m ~= n
    error('Input matrix must be square!');
end
tol = eps/4;

A = full(A);

tstart1 = tic;
[H,Q] = HW8_HHrdcUH(A);
tend1 = toc(tstart1);
fprintf('Phase I: matrix reduced to a similar upper Hessenberg in %.2f secs.\n',tend1);

maxiter = n*4;
q = 0;
tmpv = randn(length(A),1);  tmpv = tmpv/norm(tmpv);
if isreal(A) && norm(A*tmpv-(tmpv'*A)')/norm(A,'fro') >= 4*eps
    maxq = n-2;     real_nonsymm = true;
else
    maxq = n-1;     real_nonsymm = false;
end
iter = 1;
tstart2 = tic;
while q < maxq && iter <= maxiter
    
    for k = 1 : n-1
        if abs(H(k+1,k)) <= tol*(abs(H(k,k))+abs(H(k+1,k+1)))
            H(k+1,k) = 0;
        end
    end
    
    oldq = q;
    for j = n-oldq : -1 : 2
        if H(j,j-1) == 0 || (j > 2 && H(j,j-1) ~= 0 && H(j-1,j-2) == 0 && real_nonsymm)
            q = q + 1;
        else
            break;
        end
    end
    subdgH1n2 = diag(H(1:n-q,1:n-q),-1);
    p = find(subdgH1n2 == 0,1,'last');
    if isempty(p),  p = 0;  end
    sizeH22 = n-p-q;
    if q < maxq
        if sizeH22 >= 2
            evs = eig(H(n-q-1:n-q,n-q-1:n-q));
        else
            evs = H(n-q,n-q);
        end
        if isreal(evs) || ~real_nonsymm
            [~,idx] = min(abs(evs-H(n-q,n-q)));
            [H(p+1:n-q,p+1:n-q),GCS] = HW8_SingleShiftedQRstep(H(p+1:n-q,p+1:n-q),evs(idx));
            for k = 1 : sizeH22-1
                Gt = [conj(GCS(1,k)) -GCS(2,k); conj(GCS(2,k)) conj(GCS(1,k))];
                Q(:,p+k:p+k+1) = Q(:,p+k:p+k+1)*Gt;
                H(1:p,p+k:p+k+1) = H(1:p,p+k:p+k+1)*Gt;
                H(p+k:p+k+1,n-q+1:end) = Gt'*H(p+k:p+k+1,n-q+1:end);
            end
        else
            [H(p+1:n-q,p+1:n-q),HVs] = HW8_DoubleShiftedQRstep(H(p+1:n-q,p+1:n-q),H(n-q-1:n-q,n-q-1:n-q));
            for k = 1 : sizeH22-2
                Q(:,p+k:p+k+2) = Q(:,p+k:p+k+2)-Q(:,p+k:p+k+2)*(2*HVs(:,k))*HVs(:,k)';
                H(1:p,p+k:p+k+2) = H(1:p,p+k:p+k+2)-H(1:p,p+k:p+k+2)*(2*HVs(:,k))*HVs(:,k)';
                H(p+k:p+k+2,n-q+1:end) = H(p+k:p+k+2,n-q+1:end)-(2*HVs(:,k))*(HVs(:,k)'*H(p+k:p+k+2,n-q+1:end));
            end
            Q(:,n-q-1:n-q) = Q(:,n-q-1:n-q)-Q(:,n-q-1:n-q)*(2*HVs(1:2,end))*HVs(1:2,end)';
            H(1:p,n-q-1:n-q) = H(1:p,n-q-1:n-q)-H(1:p,n-q-1:n-q)*(2*HVs(1:2,end))*HVs(1:2,end)';
            H(n-q-1:n-q,n-q+1:end) = H(n-q-1:n-q,n-q+1:end)-(2*HVs(1:2,end))*(HVs(1:2,end)'*H(n-q-1:n-q,n-q+1:end));
        end
    end
    if iter == 1 || mod(iter,100) == 0 || q >= maxq
        fprintf('Iteration %d: %d eigenvalues have converged.\n',iter,q);
    end
    if q >= maxq
        fprintf('All %d eigenvalues have converged in Schur form in %d iterations.\n',n,iter);
        break;
    end
    iter = iter + 1;
end

tend2 = toc(tstart2);

if q < maxq
    fprintf('Maximum iteration (%d) reached; %d eigenvalues have converged.\n',maxiter,q);
end

fprintf('Phase II: the shifted QR iteration took %.2f seconds.\n',tend2);
end

function [H,V] = HW8_simreductHess(A)

[m,n] = size(A);
if m ~= n
    error('A must be a square matrix.');
end
tmpv = randn(length(A),1);  tmpv = tmpv/norm(tmpv);
if norm(A*tmpv-(tmpv'*A)')/norm(A,'fro') <= 2*eps
    symm = true;
else
    symm = false;
end
V = zeros(n,n-2);
for k = 1 : n-2
    xk = A(k+1:n,k);
    signxk1 = sign(xk(1));
    if signxk1 == 0,    signxk1 = 1;    end
    vk = -signxk1*norm(xk)*eye(n-k,1)-xk;
    vk = vk/norm(vk);   vk = vk/norm(vk);
    A(k+1:n,k) = -signxk1*norm(xk)*eye(n-k,1);
    A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - 2*vk*(vk'*A(k+1:n,k+1:n));
    A(:,k+1:n) = A(:,k+1:n) - 2*(A(:,k+1:n)*vk)*vk';
    V(k+1:n,k) = vk;
end
if ~symm,   H = A;
else,       H = tril(A,1);
end

end

function [H,Q] = HW8_SingleShiftedQRstep(H0,evs)
    [m,n] = size(H0);
    mu = evs*eye(m);
    Q=eye(m);
    R0=(H0-mu);
    GCS=zeros(n-1,2);
    for k=1:n-1
        [c,s] = givens(R0(k:k+1,k));
        G=[c,-s;s,c];
        R0(k:k+1,k:end)=G*R0(k:k+1,k:end);
        R0(k+1,k)=0;

        Q(1:k+1,k:k+1)=Q(1:k+1,k:k+1)*G';
        GCS(k,:)=[c,s];
    end
    H=R0*Q+mu;
    for j=1:n
        for i=j:n
            if abs(H(i,j))<10^-16
                H(i,j)=0;
            end
        end
    end
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

function [H,Q] = HW8_HHrdcUH(A)

% Householder similar reduction of a square matrix A into upper Hessenberg
% 
% Input:  A is the input square matrix, real or complex
% Output: H is an upper Hessenberg, similar to A
%         Q is orthogonal/unitary, such that Q'*A*Q = H numerically
% 
% Copyright (c) 2017, F. Xue
%

[m,n] = size(A);
if m ~= n
    error('A must be an mxn matrix where n >= n.');
end
V = zeros(n,n-2);
for k = 1 : n-2
    xk = A(k+1:n,k);
    signxk1 = sign(xk(1));
    if signxk1 == 0,    signxk1 = 1;    end
    vk = signxk1*norm(xk)*eye(n-k,1)+xk;
    vk = vk/norm(vk);   vk = vk/norm(vk);
    A(k+1:n,k) = -signxk1*norm(xk)*eye(n-k,1);
    A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - 2*vk*(vk'*A(k+1:n,k+1:n));
    A(:,k+1:n) = A(:,k+1:n) - 2*(A(:,k+1:n)*vk)*vk';
    V(k+1:n,k) = vk;
end

H = A;

Q = eye(n,n);
for k = n-2:-1:1
    Q(k+1:n,k+1:n) = Q(k+1:n,k+1:n)-2*V(k+1:n,k)*(V(k+1:n,k)'*Q(k+1:n,k+1:n));
end

end

function [Hnew,HVs] = HW8_DoubleShiftedQRstep(H,Mu2)

[m2,n2] = size(Mu2);
if m2 ~= 2 || n2 ~= 2 || ~isreal(Mu2)
    error('Wrong input Mu2: it must be a 2x2 real matrix.');
end
[m,n] = size(H);
if m ~= n
    error('Input H is not a square matrix!');
end
if nnz(tril(H,-2)) > 0
    error('Input H is not upper Hessenberg!');
end
if ~isreal(H)
    error('H is not real!');
end

HVs = zeros(3,n-1);

s = Mu2(1,1)+Mu2(2,2);
t = Mu2(1,1)*Mu2(2,2)-Mu2(1,2)*Mu2(2,1);

x = H(1,1)^2+H(1,2)*H(2,1)-s*H(1,1)+t;
y = H(2,1)*(H(1,1)+H(2,2)-s);
z = H(2,1)*H(3,2);
for k = 0 : n-3
    if y == 0 && z == 0
        continue;
    end
    signx = sign(x);
    if signx == 0,   signx = 1;  end
    vk = signx*norm([x; y; z])*eye(3,1) + [x; y; z];
    vk = vk/norm(vk);   vk = vk/norm(vk);
    HVs(:,k+1) = vk;
    q = max([1 k]);
    H(k+1:k+3,q:n) = H(k+1:k+3,q:n) - (2*vk)*(vk'*H(k+1:k+3,q:n));
    r = min([k+4 n]);
    H(1:r,k+1:k+3) = H(1:r,k+1:k+3) - (H(1:r,k+1:k+3)*(2*vk))*vk';
    x = H(k+2,k+1);
    y = H(k+3,k+1);
    if k < n-3
        z = H(k+4,k+1);
    end
end
if y ~= 0
    signx = sign(x);
    if signx == 0,  signx = 1;  end
    vk = signx*norm([x; y])*eye(2,1) + [x; y];
    vk = vk/norm(vk);   vk = vk/norm(vk);
    HVs(1:2,n-1) = vk;
    H(n-1:n,n-2:n) = H(n-1:n,n-2:n)-(2*vk)*(vk'*H(n-1:n,n-2:n));
    H(1:n,n-1:n) = H(1:n,n-1:n)-(H(1:n,n-1:n)*(2*vk))*vk';
end
Hnew = triu(H,-1);
end