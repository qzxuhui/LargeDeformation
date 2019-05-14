clc
clear
E=eye(3,3);
syms J
syms lambda mu
syms s1 s2 s3 s4 s5 s6
syms F11 F12 F13
syms F21 F22 F23
syms F31 F32 F33

%% bmat
F=[F11 F12 F13
   F21 F22 F23
   F31 F32 F33];
Ft=[F11 F21 F31
    F12 F22 F32
    F13 F23 F33];
B=F*Ft

%% ddsdde
S=[
    s1 s4 s5
    s4 s2 s6
    s5 s6 s3
];
for i=1:1:3
    for j=1:1:3
        for k=1:1:3
            for l=1:1:3
                a=transfer_index_to_vogit(i,j);
                b=transfer_index_to_vogit(k,l);
                if a==0 || b==0
                    continue;
                end
                t1=lambda*E(i,j)*E(k,l)+mu*(E(i,k)*E(j,l)+E(i,l)*E(k,j));
                t1=t1/J;
                t2=1/2*(S(i,k)*E(j,l)+S(j,l)*E(i,k)+S(i,l)*E(j,k)+S(j,k)*E(i,l));
                DDSDDE(a,b)=t1+t2;
            end
        end
    end
end

DDSDDE

function a=transfer_index_to_vogit(i,j)
a=0;
if i==1 && j==1
    a=1;
end
if i==2 && j==2
    a=2;
end
if i==3 && j==3
    a=3;
end
if i==1 && j==2
    a=4;
end
if i==1 && j==3
    a=5;
end
if i==2 && j==3
    a=6;
end
end