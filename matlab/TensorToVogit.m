clc
clear
close all
syms lambda mu p q
E=eye(3);
dim=3;
for i=1:1:dim
    for j=1:1:dim
        for l=1:1:dim
            for k=1:1:dim
                a=0;
                b=0;
                a=transfer_vogit(i,j);
                b=transfer_vogit(k,l);
                if (a ==0 || b==0)
                    continue;
                end
                C(a,b)=2*p*E(i,k)*E(j,l)+q*E(i,j)*E(k,l);
            end
        end
    end
end
C
function a=transfer_vogit(i,j)
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
if i==2 && j==3
    a=4;
end
if i==1 && j==3
    a=5;
end
if i==1 && j==2
    a=6;
end

end