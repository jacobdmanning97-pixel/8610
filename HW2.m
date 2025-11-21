%problem 1
list=[];

n=9;
list=[list lisn(n)];

n=19;
list=[list lisn(n)];

n=29;
list=[list lisn(n)];

n=39;
list=[list lisn(n)];

%parts a-d
disp(list(1:4,:))

%part e
disp(list(5,:))

%problem 3
for j = 5:9
    i = 0;
    n = 2^j;
    
    for k = 1:1000
        Ap=randn(n,n)/sqrt(n);
        ak = luFactor(Ap);
        am = max(abs(Ap(:)));
        p = ak/am;
        
        if p > sqrt(n)
            i=i+1;
        end
    end

    disp(i/k)

end

function l=lisn(n)
    l=zeros(5,2);
    u=linspace(-1,1,n+1);
    A=vander(u);
    x=repelem(1,n+1)';
    b=A*x;
    [Q,R]=qr(A);
    
    x1=A\b;
    x2=R\(Q'*b);
    x3=cram(b,A);
    x4=inv(A)*b;
    x5=GEnoP(A, b);
    
    [fe1,be1]=err(x1,x,b,A);
    [fe2,be2]=err(x2,x,b,A);
    [fe3,be3]=err(x3,x,b,A);
    [fe4,be4]=err(x4,x,b,A);
    [fe5,be5]=err(x5,x,b,A);
    
    l(1,:)=[fe1,be1];
    l(2,:)=[fe2,be2];
    l(3,:)=[fe3,be3];
    l(4,:)=[fe4,be4];
    l(5,:)=[fe5,be5];
end

function [fe, be]=err(xh,x,b,A)
    fe=norm(xh-x)/norm(x);
    be=norm(b-A*xh)/(norm(A)*norm(xh));
end

function X=cram(b,A)
    % Determinant of coefficient matrix
    det_A = det(A);
    
    % Solution vector
    X = zeros(size(b));
    
    % Cramer's rule
    for i = 1:size(A, 1)
        % Replace the ith column of A with B
        A_i = A;
        A_i(:, i) = b;
        
        % Calculate the determinant of A_i
        det_A_i = det(A_i);
        
        % Calculate the ith element of X
        X(i) = det_A_i / det_A;
    end
end

function x = GEnoP(A, b)
    B = [A, b];
    [n, m] = size(B);
    % Start with Forward sub
    for i = 1:n
        B(i, :) = B(i, :) ./ B(i, i);
        for k = i+1:n
            B(k, :) = (-B(i, :) * B(k, i)) + B(k, :);
        end
    end
    % Back Substitution
    for j = n-1:-1:1
        for z = j+1:1:n
            B(j, :) = (-B(z, :) * B(j, z)) + B(j, :);
        end
    end
    x = B(:, end);
end

function ak = luFactor(A)
    % Get the size of A
    [n, ~] = size(A);

    % Initialize L, U, and P
    L = eye(n);
    U = A;
    P = eye(n);
    ak = 0;

    for k = 1:n-1
        % Find the maximum element in the current column
        [~, i] = max(abs(U(k:n, k)));
        if max(abs(U(k:n, k))) > ak
            ak = max(abs(U(k:n, k)));
        end
        i = i + k - 1;

        % Swap rows i and k in U, and record the permutation in P
        U([k, i], :) = U([i, k], :);
        P([k, i], :) = P([i, k], :);
    end
end
