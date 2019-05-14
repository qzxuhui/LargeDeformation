clc
clear
close all
%% init
syms a11 a12 a21 a22
syms sigma11 sigma12 sigma21 sigma22
syms J
syms lambda mu p
A=[a11,a12;a21,a22];
Sigma=[sigma11,0;0,0];
E=eye(2);
%% C_tau
for i=1:1:2
    for j=1:1:2
        for k=1:1:2
            for l=1:1:2
                C_tau(i,j,k,l)=lambda*E(i,j)*E(k,l)+mu*(E(i,k)*E(j,l)+E(i,l)*E(k,j));
            end
        end
    end
end
%% C_sigmaT
C_sigmaT=C_tau/J;
%% pPpF
for i=1:1:2
    for j=1:1:2
        for k=1:1:2
            for l=1:1:2
                pP_pF(i,j,k,l)=C_sigmaT(i,j,k,l)+Sigma(i,l)*E(j,k);
            end
        end
    end
end
%% get enerey
sum=0;
for i=1:1:2
    for j=1:1:2
        for k=1:1:2
            for l=1:1:2
                sum=sum+1/2*(A(i,j)+A(j,i))*pP_pF(i,j,k,l)*(A(k,l));
            end
        end
    end
end
sum=expand(sum);
%% get
for i=1:1:2
    for j=1:1:2
        for k=1:1:2
            for l=1:1:2
                a=transfer_index(i,j);
                b=transfer_index(k,l);
                pP_pF_mod(i,j,k,l)=diff(diff(sum,A(i,j)),A(k,l));
            end
        end
    end
end
C_tau_mat=Transfer_Tensor_to_Matrix(C_tau)
pP_pF_mat=Transfer_Tensor_to_Matrix(pP_pF)
pP_pF_mod_mat=Transfer_Tensor_to_Matrix(pP_pF_mod)/2

%% function
function  mat=Transfer_Tensor_to_Matrix(ten)

for i=1:1:2
    for j=1:1:2
        for k=1:1:2
            for l=1:1:2
                a=transfer_index(i,j);
                b=transfer_index(k,l);
                mat(a,b)=ten(i,j,k,l);
            end
        end
    end
end

end

function a=transfer_index(i,j)
% a=0;
if i==1 && j==1
    a=1;
end
if i==2 && j==2
    a=2;
end
if i==1 && j==2
    a=3;
end
if i==2 && j==1
    a=4;
end

end