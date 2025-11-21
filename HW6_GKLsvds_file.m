disp("PicA")
picA = double(imread('picA.jpg'));
[pic, pich] = comp(picA,160);
figure(1); imshow(uint8(pic)); axis equal;
figure(2); imshow(uint8(pich)); axis equal;

[~, pich] = comp(picA,320);
figure(3); imshow(uint8(pich)); axis equal;

disp("PicB")
picB = double(imread('picB.jpg'));
[pic, pich] = comp(picB,160);
figure(4); imshow(uint8(pic)); axis equal;
figure(5); imshow(uint8(pich)); axis equal;

[~, pich] = comp(picB,320);
figure(6); imshow(uint8(pich)); axis equal;

disp("PicC")
picC = double(imread('picC.jpg'));
[pic, pich] = comp(picC,160);
figure(7); imshow(uint8(pic)); axis equal;
figure(8); imshow(uint8(pich)); axis equal;

[~, pich] = comp(picC,320);
figure(9); imshow(uint8(pich)); axis equal;

disp("PicD")
picD = double(imread('picD.jpg'));
[pic, pich] = comp(picD,160);
figure(10); imshow(uint8(pic)); axis equal;
figure(11); imshow(uint8(pich)); axis equal;

[~, pich] = comp(picD,320);
figure(12); imshow(uint8(pich)); axis equal;

disp("PicE")
picE = double(imread('picE.jpg'));
[pic, pich] = comp(picE,160);
figure(13); imshow(uint8(pic)); axis equal;
figure(14); imshow(uint8(pich)); axis equal;

[~, pich] = comp(picE,320);
figure(15); imshow(uint8(pich)); axis equal;

disp("PicF")
picF = double(imread('picF.jpg'));
[pic, pich] = comp(picF,160);
figure(16); imshow(uint8(pic)); axis equal;
figure(17); imshow(uint8(pich)); axis equal;

[~, pich] = comp(picF,320);
figure(18); imshow(uint8(pich)); axis equal;

disp("PicG")
picG = double(imread('picG.jpg'));
[pic, pich] = comp(picG,160);
figure(19); imshow(uint8(pic)); axis equal;
figure(20); imshow(uint8(pich)); axis equal;

[~, pich] = comp(picG,320);
figure(21); imshow(uint8(pich)); axis equal;

 function [pic, pich] = comp(pic,rk)
     tic; [Us1,Ss1,Vs1] = HW6_GKLsvds(pic(:,:,1),rk); toc;
     tic; [Us2,Ss2,Vs2] = HW6_GKLsvds(pic(:,:,2),rk); toc;
     tic; [Us3,Ss3,Vs3] = HW6_GKLsvds(pic(:,:,3),rk); toc;
     tic; [U1,S1,V1] = svd(pic(:,:,1),0); toc;
     tic; [U2,S2,V2] = svd(pic(:,:,2),0); toc;
     tic; [U3,S3,V3] = svd(pic(:,:,3),0); toc;
    
    whos pic
    whos Us1 Vs1 Us2 Vs2 Us3 Vs3
    whos U1 V1 U2 V2 U3 V3
    
     pich = zeros(size(pic));
     pich(:,:,1) = Us1*Ss1*Vs1';
     pich(:,:,2) = Us2*Ss2*Vs2';
     pich(:,:,3) = Us3*Ss3*Vs3';
    
     disp([norm(pich(:,:,1)-pic(:,:,1),'fro')/norm(pic(:,:,1),'fro') ...
     norm(pich(:,:,2)-pic(:,:,2),'fro')/norm(pic(:,:,2),'fro') ...
     norm(pich(:,:,3)-pic(:,:,3),'fro')/norm(pic(:,:,3),'fro')]);
 end

function [U,S,V] = HW6_GKLsvds(A,k)

% The Golub-Kahan-Lanczos bidiagonalization
%
% Input:
% A     The matrix for which we are computing largest singular values
%       should be large to show the competitiveness of iterative method
% k     The number of singular values wanted
%
% Output:
% S     The k by k diagonal matrix of approximate dominant singular values
% U,V   The k approximate left and right singular vectors
%
% Copyright (c) F. Xue 10/21/2017

[m,n] = size(A);
v_k = randn(n,1);
v_k = v_k/norm(v_k);
maxiter = min([m n max([ceil(1.2*k) k+5])]);
beta_km1 = 0;
u_km1 = zeros(m,1);
alpha_all = zeros(maxiter,1);
beta_all = zeros(maxiter+1,1);
U = zeros(m,maxiter);
V = zeros(n,maxiter+1);
V(:,1) = v_k;

for iter = 1 : maxiter
    u_k = A*v_k-beta_km1*u_km1;
    
    for jj = 1 : iter-1
        u_k = u_k - U(:,jj)*(U(:,jj)'*u_k);
    end
    %u_k = u_k - U(:,1:iter-1)*(U(:,1:iter-1)'*u_k);
    
    alpha_k = norm(u_k);
    u_k = u_k/alpha_k;
    v_kp1 = (u_k'*A)'-alpha_k*v_k;
    
    for jj = 1 : iter
        v_kp1 = v_kp1 - V(:,jj)*(V(:,jj)'*v_kp1);
    end
    %v_kp1 = v_kp1 - V(:,1:iter)*(V(:,1:iter)'*v_kp1);
    
    beta_k = norm(v_kp1);
    v_kp1 = v_kp1/beta_k;
    
    alpha_all(iter) = alpha_k;
    beta_all(iter+1) = beta_k;
    
    U(:,iter) = u_k;
    V(:,iter+1) = v_kp1;
    
    u_km1 = u_k;
    v_k = v_kp1;
    beta_km1 = beta_k;
end

B = spdiags([alpha_all beta_all(1:end-1)],0:1,maxiter,maxiter);
[Us,S,Vs] = svd(full(B));
U = U*Us(:,1:k);
V = V(:,1:maxiter)*Vs(:,1:k);
S = S(1:k,1:k);

end
