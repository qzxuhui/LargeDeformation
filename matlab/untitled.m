clc
clear
close all
lambda_2=linspace(1,2,300);
eps=0.0001;
A=lambda_2.^2;
J=exp(eps*(1-A));
sigma_c=1./(J.^2).*(1./A.*J-A);
sigma_u=(1./A-A);
figure
hold on
plot(lambda_2,sigma_u);
plot(lambda_2,sigma_c);
ylabel('$$\frac{\sigma_{11}}{\mu_0}$$','interpreter','latex');
xlabel('$$\lambda_2$$','interpreter','latex');

dsigma_c_dlambda=Diff(sigma_c,lambda_2);
dsigma_u_dlambda=Diff(sigma_u,lambda_2);
plot(lambda_2,dsigma_c_dlambda)
plot(lambda_2,dsigma_u_dlambda)

figure
plot(lambda_2,J)
function [dydx]=Diff(y,x)
n=length(y)
dydx=zeros(n,1);
for i=1:1:n-1
    dy=y(i+1)-y(i);
    dx=x(i+1)-x(i);
    dydx(i)=dy/dx;
end

end