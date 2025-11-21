A1 = randn(2^20,15);
u = (-1:2/40:1)';
A2 = u.^(0:23);
A3 = u.^(0:40);

%problem 2
l1 = list(A1);

l1 = ["","A1",""; l1];

l2 = list(A2);

l2 = ["","A2",""; l2];

l3 = list(A3);

l3 = ["","A3",""; l3];

%problem 4
[nor1, sol1] = LS(A2);
[nor2, sol2] = LS(A3);

l4=["","Calculated", "MatLab"; "A2", nor1, sol1; "A3", nor2, sol2];

l1
l2
l3
l4

function l = list(A)
    l = zeros(4,3);
    
    [Q, R] = CGS(A);
    [n1, n2] = normp2(A, Q, R);
    l(1,:) = [1, n1, n2];
    
    [Q, R] = MGS(A);
    [n1, n2] = normp2(A, Q, R);
    l(2,:) = [2, n1, n2];
    
    [Q, R] = MGSwRO(A);
    [n1, n2] = normp2(A, Q, R);
    l(3,:) = [3, n1, n2];
    
    [Q, R, V] = HQR(A);
    [n1, n2] = normp2(A, Q, R);
    l(4,:) = [4, n1, n2];
end 

function [n1, n2] = normp2(A, Q, R)
    [m, n] = size(Q);
    n1 = norm(A-Q*R,"fro")/norm(A,"fro");
    n2 = norm(Q'*Q - eye(n));
end

function [nor, sol] = LS(A)
    [Q, R, V] = HQR(A);
    [m, n] = size(A);
    b = -cumprod(-ones(1,m))';
    bp = -cumprod(-ones(1,m))';
    V = fliplr(V);
    
    for i = 1:n
        b = b - 2*V(:,i)*(V(:,i)'*b);
    end
    
    b = b(1:n,:);
    
    xh = R\b;
    xt = A\bp;

    nor = norm(bp - A*xh)/(norm(A)*norm(xh));
    sol = norm(bp - A*xt)/(norm(A)*norm(xt));
end

function [Q, R] = CGS(A)
    [n, m] = size(A); 
    Q = zeros(n, m);
    R = zeros(m, m);
    
    R(1, 1) = norm(A(:, 1));
    Q(:, 1) = A(:, 1) / R(1, 1);
    
    for j = 2:m
        Q(:, j) = A(:, j);
        for i = 1:j-1
            Q(:, j) = Q(:, j) - Q(:, i) * (Q(:, i)' * A(:, j));
            R(i, j) = Q(:, i)' * A(:, j);
        end
        R(j, j) = norm(Q(:, j));
        Q(:, j) = Q(:, j) / R(j, j);
    end
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

% function [Q, R] = MGS(A)
%     [n, m] = size(A); 
%     Q = zeros(n, m);
%     R = zeros(m, m);
% 
%     R(1, 1) = norm(A(:, 1));
%     Q(:, 1) = A(:, 1) / R(1, 1);
% 
%     for j = 2:m
%         Q(:, j) = A(:, j);
%         for i = 1:j-1
%             Q(:, j) = Q(:, j) - Q(:, i) * (Q(:, i)' * A(:, j));
%             R(i, j) = Q(:, i)' * A(:, j);
%         end
%         R(j, j) = norm(Q(:, j));
%         Q(:, j) = Q(:, j) / R(j, j);
%     end
% end

function [Q, R] = MGSwRO(A)
    [n, m] = size(A); 
    Q = zeros(n, m);
    R = zeros(m, m);
    
    R(1, 1) = norm(A(:, 1));
    Q(:, 1) = A(:, 1) / R(1, 1);
    
    for j = 2:m
        Q(:, j) = A(:, j);
        for i = 1:j-1
            Q(:, j) = Q(:, j) - Q(:, i) * (Q(:, i)' * A(:, j));
            R(i, j) = Q(:, i)' * A(:, j);
        end
        for i = 1:j-1
            Q(:, j) = Q(:, j) - Q(:, i) * (Q(:, i)' * A(:, j));
            R(i, j) = R(i, j) + Q(:, i)' * A(:, j);
        end
        R(j, j) = norm(Q(:, j));
        Q(:, j) = Q(:, j) / R(j, j);
    end
end

function [Q, R, V] = HQR(A)
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
