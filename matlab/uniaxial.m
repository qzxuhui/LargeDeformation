function y=uniaxial(x,mu_0,lambda_0)
% variable need to solve
lambda_1=x(1);
lambda_2=x(2);

% jacobian
J=lambda_1*lambda_2;

% properties
lambda=lambda_0;
mu=mu_0-lambda*log(J);

% stress
stress_11=mu_0/J*(lambda_1*lambda_1)-mu/J;
stress_22=mu_0/J*(lambda_2*lambda_2)-mu/J;

% equations
y=zeros(2);
y(2)=stress_11/2+mu/J;
y(3)=stress_22;
end