format long
rng('default');
H0 = triu(randn(7,7),-1);

flag = true;
m = zeros(8,2);
[l,n]=normf(H0,3,flag);
m(1,:)=[l,n];
[l,n]=normf(H0,30,flag);
m(2,:)=[l,n];
[l,n]=normf(H0,300,flag);
m(3,:)=[l,n];
[l,n]=normf(H0,3000,flag);
m(4,:)=[l,n];

flag = false;
[l,n]=normf(H0,3,flag);
m(5,:)=[l,n];
[l,n]=normf(H0,30,flag);
m(6,:)=[l,n];
[l,n]=normf(H0,300,flag);
m(7,:)=[l,n];
[l,n]=normf(H0,3000,flag);
m(8,:)=[l,n];

m

[V,d]=eig(H0,'vector');
d=sort(unique(abs(d)),'descend');

l=zeros(size(d)-1);

for i = 1:size(d)-1
    l(i)=d(i+1)/d(i);
end
max(l)

function [l,n] = normf(A,k,flag)
    H0=SI(A,k,flag);
    H1=UQR(A,k,flag);
    l=norm(H0-H1,'fro')/norm(H0,'fro');
    
    [V0,d0]=eig(H0,'vector');
    [V1,d1]=eig(H1,'vector');
    d0=sort(d0);
    d1=sort(d1);
    n=norm(d0-d1)/norm(d0);
end

function H = SI(A,k,flag)
    Q=eye(size(A));
    for i=1:k
        Z=A*Q;
        [Q,R]=qr(Z);
        if flag == true
            [Q,R]=PP(Q,R);
        end
    end
    H=Q'*A*Q;
end

function H = UQR(A,k,flag)
    H=A;
    for i=1:k
        [Q,R]=qr(H);
        if flag == true
            [Q,R]=PP(Q,R);
        end
        H=R*Q;
    end
end

function [Q,R]=PP(Q, R)
    for i=1:size(R)
        if R(i,i)<0
            R(i,:)=-R(i,:);
            Q(:,i)=-Q(:,i);
        end
    end
end