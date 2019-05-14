% progressive analysis to found when instabilty occurs

clc
clear
close all
%% properties
mu_0=2;
lambda_0=100000;
lambda_2=linspace(1,10,100);

%% solver
% options = optimoptions('fsolve','algorithm','Levenberg-Marquardt',...
%     'Display','iter','MaxFunctionEvaluations',10000,'MaxIterations',10000);
% solution_0=[0.5,2];
% solve=@(x)uniaxial(x,mu_0,lambda_0);
% sol=fsolve(solve,solution_0,options)
% % why residual is 3?

%% compressible analysis result
clambda_2=lambda_2;
clambda_1=1./clambda_2.*exp(mu_0/lambda_0*(1-clambda_2.^2));
cJ=clambda_1.*clambda_2;
cmu=mu_0-lambda_0.*log(cJ);
csigma_1=mu_0./cJ.*clambda_1.^2-cmu./cJ;
csigma_2=mu_0./cJ.*clambda_2.^2-cmu./cJ;
%% uncompressible analysis result
ulambda_2=linspace(1,10,100);
ulambda_1=1./ulambda_2;
uJ=ulambda_1.*ulambda_2;
umu=mu_0-lambda_0.*log(uJ);
up=mu_0*lambda_2.^2;
usigma_1=mu_0.*ulambda_1.^2-up;
usigma_2=mu_0.*ulambda_2.^2-up;
%% compressible plot
figure
set(gcf,'position',[50,50,1500,800]);
subplot(2,2,1)
plot(clambda_1,clambda_2,'linewidth',2)
xlabel('$$\lambda_1$$','interpreter','latex');
ylabel('$$\lambda_2$$','interpreter','latex');
set(gca,'fontsize',18)
subplot(2,2,2)
plot(clambda_1,cJ,'linewidth',2)
xlabel('$$\lambda_1$$','interpreter','latex');
ylabel('$$J$$','interpreter','latex');
set(gca,'fontsize',18)
subplot(2,2,3)
hold on
plot(clambda_1,csigma_1,'linewidth',2)
plot(clambda_1,-2*cmu./cJ,'linewidth',2)
xlabel('$$\lambda_1$$','interpreter','latex');
ylabel('$$\sigma_1$$','interpreter','latex');
set(gca,'fontsize',18)
subplot(2,2,4)
plot(clambda_1,csigma_2,'linewidth',2)
xlabel('$$\lambda_1$$','interpreter','latex');
ylabel('$$\sigma_2$$','interpreter','latex');
set(gca,'fontsize',18)
%% uncompressible plot
figure
set(gcf,'position',[50,50,1500,800]);
subplot(2,2,1)
plot(ulambda_1,ulambda_2,'linewidth',2)
xlabel('$$\lambda_1$$','interpreter','latex');
ylabel('$$\lambda_2$$','interpreter','latex');
set(gca,'fontsize',18)
subplot(2,2,2)
plot(ulambda_1,uJ,'linewidth',2)
xlabel('$$\lambda_1$$','interpreter','latex');
ylabel('$$J$$','interpreter','latex');
set(gca,'fontsize',18)
subplot(2,2,3)
hold on
plot(ulambda_1,usigma_1,'linewidth',2)
plot(ulambda_1,-2*umu./uJ,'linewidth',2)
xlabel('$$\lambda_1$$','interpreter','latex');
ylabel('$$\sigma_1$$','interpreter','latex');
set(gca,'fontsize',18)
subplot(2,2,4)
plot(ulambda_1,usigma_2,'linewidth',2)
xlabel('$$\lambda_1$$','interpreter','latex');
ylabel('$$\sigma_2$$','interpreter','latex');
set(gca,'fontsize',18)
close all
%% compare
figure
title('x')
set(gcf,'position',[50,50,1500,800]);

subplot(2,3,1)
hold on
plot(lambda_2,clambda_1,'s','linewidth',2);
plot(lambda_2,ulambda_1,'.','linewidth',2);
xlabel('$$\lambda_2$$','interpreter','latex');
ylabel('$$\lambda_1$$','interpreter','latex');
set(gca,'fontsize',18)

subplot(2,3,2)
hold on
plot(lambda_2,csigma_1,'s','linewidth',2);
plot(lambda_2,usigma_1,'.','linewidth',2);
xlabel('$$\lambda_2$$','interpreter','latex');
ylabel('$$\sigma_1$$','interpreter','latex');
set(gca,'fontsize',18)

subplot(2,3,3)
hold on
plot(lambda_2,csigma_2,'s','linewidth',2);
plot(lambda_2,usigma_2,'.','linewidth',2);
xlabel('$$\lambda_2$$','interpreter','latex');
ylabel('$$\sigma_2$$','interpreter','latex');
set(gca,'fontsize',18)

subplot(2,3,4)
hold on
plot(lambda_2,cmu,'s','linewidth',2);
plot(lambda_2,umu,'.','linewidth',2);
xlabel('$$\lambda_2$$','interpreter','latex');
ylabel('$$\mu$$','interpreter','latex');
set(gca,'fontsize',18)

subplot(2,3,5)
hold on
plot(lambda_2,cJ,'s','linewidth',2);
plot(lambda_2,uJ,'.','linewidth',2);
xlabel('$$\lambda_2$$','interpreter','latex');
ylabel('$$J$$','interpreter','latex');
set(gca,'fontsize',18)

subplot(2,3,6)
hold on
plot(lambda_2,cmu./cJ,'s','linewidth',2);
plot(lambda_2,umu./uJ,'.','linewidth',2);
xlabel('$$\lambda_2$$','interpreter','latex');
ylabel('$$\frac{\mu}{J}$$','interpreter','latex');
set(gca,'fontsize',18)
%%
figure
hold on
plot(clambda_1,csigma_1/2+cmu./cJ);
plot(ulambda_1,usigma_1/2+umu./uJ);
plot(ulambda_1,0*ulambda_1,'k--');