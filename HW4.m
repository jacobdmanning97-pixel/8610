format long

u = (-1:2/40:1)';
A2 = u.^(0:23);
A3 = u.^(0:40);
U = randn(1024,10);
A4 = U*randn(10,15);

%Problem 11.3
t = linspace(0,1,50);
A1 = fliplr(vander(t));
A1 = A1(:,1:12);
b = cos(4*t)';

x1 = (A1'*A1)\A1'*b;

[Q,R] = MGS(A1);
x2 = R\Q'*b;

[Q,R] = HQR(A1);
x3 = R\Q'*b;

[Q,R] = qr(A1);
x4 = R\Q'*b;

x5 = A1\b;

[U,S,V] = svd(A1,"econ");
x6 = V*(S\U')*b;

x = [x1 x2 x3 x4 x5 x6];

%problem 3
n2 = prob3(A2);
n3 = prob3(A3);
n4 = prob3(A4);

n = [n2 n3 n4];

%problem 4
[m, n] = size(A2);
B2 = [zeros(n,n); A2];

[m, n] = size(A3);
B3 = [zeros(n,n); A3];

[Q2, R2, P2] = qr(B2, 0);
[Q3, R3, P3] = qr(B3, 0);

function n = prob3(A)
    [Q, R, P] = HQRwP(A);
    [Qm, Rm, Pm] = qr(A, 0);
    n = norm(R-Rm);
end

function [Q, R] = MGS(A)
    [n, m] = size(A); 
    Q = A;
    R = zeros(m, m);
    
    for k = 1:m
        for i = 1:k-1
            R(i,k) = Q(:,i)'*Q(:,k);
            Q(:,k) = Q(:,k) - R(i,k)*Q(:,i);
        end
        R(k,k) = norm(Q(:,k));
        Q(:,k) = Q(:,k)/R(k,k);
    end
end

function [Q, R] = HQR(A)
    [m,n] = size(A);
    R = A;
    V = zeros(m,n);
    Q = [eye(n);zeros(m-n,n)];

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
        V(:,k) = [zeros(k-1,1);v];
    end

    R = R(1:n,1:n);
    V = fliplr(V);

    for i = 1:n
        Q = Q - 2*V(:,i)*(V(:,i)'*Q);
    end

end

function [Q, R, P] = HQRwP(A)
    [m,n] = size(A);
    R = A;
    P = eye(n);
    V = zeros(m,n);
    Q = [eye(n);zeros(m-n,n)];

    for k = 1:n
        PP = eye(n);
        y = zeros(m-k,1);
        for j = k:n
            y(j,1) = norm(R(k:m,j));
        end

        [B,I] = sort(y,'descend');

        z = R(:,k);
        R(:,k) = R(:,I(1));
        R(:,I(1)) = z;

        z = PP(:,k);
        PP(:,k) = PP(:,I(1));
        PP(:,I(1)) = z;
        P = P*PP;

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
        V(:,k) = [zeros(k-1,1);v];
    end

    R = R(1:n,1:n);
    V = fliplr(V);

    for i = 1:n
        Q = Q - 2*V(:,i)*(V(:,i)'*Q);
    end

end
