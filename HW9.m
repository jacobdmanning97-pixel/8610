rng('default'); n = 1024; A = randn(n,n); [A,R] = qr(A);
 Ahat = A+1.2*eye(n); E = randn(n,n); E = E+E';
 B = (A+A')/2; B = B+1e-4*E; Bhat = B+1.01*eye(n);

 D=eig(B);
 D=sort(D);

 F=eig(Bhat);
 F=sort(F);

 ca=cond(A);
 cahat=cond(Ahat);
 cb=cond(B);
 cbhat=cond(Bhat);

 L=eig(A);
 M=eig(Ahat);

 hold on
 plot(real(M),imag(M),'o')
 plot(real(L),imag(L),'o')
 plot(real(D),imag(D),'o')
 plot(real(F),imag(F),'o')
 hold off

  f = ones(n,1); m = n-1; restart = 1; tol = 1e-12;
 [x1,flag1,relres1,iter1,resvec1] = gmres(A,f,m,tol,1);
 semilogy(resvec1/norm(f),'ro'); hold on;
 [x2,flag2,relres2,iter2,resvec2] = gmres(Ahat,f,m,tol,1);
 semilogy(resvec2/norm(f),'go'); hold on;
 [x3,flag3,relres3,iter3,resvec3] = minres(B,f,tol,m);
 semilogy(resvec3/norm(f),'bo'); hold on;
 [x4,flag4,relres4,iter4,resvec4] = pcg(Bhat,f,tol,m);
 semilogy(resvec4/norm(f),'ko'); hold on;
 legend('A by gmres','Ahat by gmres','B by minres','Bhat by cg');

