function [Omega, Damping, VibrationMode]=GenerateMode(A, B, G, T)

%% input
%      A: State Matrix, 2n*2n. n is DOFs, 2*n is the number of state variables
%      B: Input Matrix, 2n*p. p is the number of input
%      G: Output Matrix, q*2n. q is the number of output
%      T: Sampling time, real number
%% output
%      Omega: Frequency Vector, n*1
%      Damping: Damping Ratio Vector, n*1
%      VibrationMode: Vibration Mode, n*n matrix

[VM_notSorted, Omega_notSorted]=eig(A);
Omega_notSorted=diag(logm(Omega_notSorted)./T);

%% Sorting frequency according to Imaginary part of Omega_notSorted
PositiveImagNumber=find(imag(Omega_notSorted)>0);
% PositiveImagNumber is the serial number of elements, no order.
[Dummy,SortNumber]=sort(imag(Omega_notSorted(PositiveImagNumber))); % Sorting frequency descending order
% "SortNumber" is a relative sequence of Omega_notSorted's element which has a positive imaginary, in descending order considering imaginary 

%% Inducing frequency in ascending order
SortResult=PositiveImagNumber(SortNumber); % SortResult is the descending sequence of all elements with positive imaginary.

for i=1:length(SortResult)
	lambda(i)=Omega_notSorted(SortResult(i));% Sorting frequency in descending order of imaginary
	VibrationMode(:,i)=VM_notSorted(:,SortResult(i));% Sorting eigenvector in descending order of imaginary
end

Omega=abs(lambda); % Frequency is module of lambda
Damping=-real(lambda)./Omega;% Damping ratio is real part of lambda divided by frequency