function [A,B,G]=Era(h,n,s,r)

%% input
%      h: Impulse response
%      n: Degree of freedoms(Only if structural system)
%      r: number of rows
%      s: number of columns
%% output
%      A: State Matrix, 2n*2n. n is DOFs, 2*n is the number of state variables
%      B: Input Matrix, 2n*p. p is the number of input
%      G: Output Matrix, q*2n. q is the number of output

[MM,NN]=size(h);

%% Generate Hankel Matrix H0 & H1
H0=cell(r,s);
H1=cell(r,s);

k=1; % Hankel Matrix H0
for i=1:r
	for j=1:s
		H0{i,j}=h(:,k+i+j-2);
	end
end

k=2; % Hankel Matrix H1
for i=1:r
	for i=1:s
		H1{i,j}=h(:,k+i+j-2);
	end
end
% Transform cell to matrix
H0=cell2mat(H0);
H1=cell2mat(H1);

%% Singular Value Decomposition
[U,S,V] = svd(H0,0);
% Extract first n*n part of matrix
U = U(:,1:2*n);
S = S(1:2*n,1:2*n);
V = V(:,1:2*n);

%% Generate Auxiliary Matrix Ey and Ep
ny=MM;np=1; % ny=MM and np=1 because ny should be equal to the number of outputs and np should be equal to the number of inputs
Ey = [eye(ny);zeros(ny*(r-1),ny)]; % Generate Ey(nyr*ny)
Ep = [eye(np);zeros(np*(s-1),np)]; % Generate Ep(nps*np)

%% Induce the A B G matrix
A = S^(-1/2) * U' *H1 * V * S^(-1/2);
B = S^(1/2) * V' * Ep;
G = Ey' * U *S^(1/2);


