%% Load input
A=load('/net/data1/riklund/flux_output/nmax30states30contour0/densities/density_0.050.dat');
sz=sqrt(size(A,1));
B=reshape(A(:,3),sz, sz);
X=reshape(A(:,1),sz,sz);
Y=reshape(A(:,2),sz,sz);

%% Perform calculations and plot
C=fft2(B);
%Kx=fft2(X);
%Ky=fft2(Y);
subplot(1,2,1);
contourf(X, Y, log10(real(C)));
title('Real part');
colorbar;0
subplot(1,2,2);
contourf(X,Y, log10(imag(C)));
title('Imag part');
colorbar;