A=load('output/nmax30states30contour0/gradients/gradient_0.100.dat');
X=reshape(A(:,1),500,500);
Y=reshape(A(:,2),500,500);
Zx=reshape(A(:,3),500,500);
Zy=reshape(A(:,4),500,500);
contourf(X,Y,log10(Zx+Zy));
colorbar;