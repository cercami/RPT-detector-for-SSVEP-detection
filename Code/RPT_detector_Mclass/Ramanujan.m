%% RAMANUJAN SPACE
function [Ramanujan_Space,divisors] = Ramanujan(P)
%P = input('Please enter a value for the periodicity!');
Ramanujan_Space = zeros(P,P);
counter = 1;
for d=1:P
    if mod(P,d)==0
        divisors(1,counter)=d;
        counter=counter+1;
    end
end

for i=1:counter-1
q=divisors(i);
    for n = 0:q-1
    k = 1:q;
    c_q = exp(1j*2*pi*k*n/q);
    temp_c_q = 0;
    for k = 1:q
        if gcd(k,q)==1
            temp_c_q = temp_c_q+c_q(k);
        end
    end
    ramanujan(n+1) = temp_c_q;
    ramanujan = real(ramanujan);
end
%% Euler totient function Phi(n)
Phi_Euler = 1;
for coprime=2:q
    f_1 = factor(coprime);
    f_2 = factor(q);
    check = ismember(f_1,f_2); 
    check = sum(check);
    if check==0
       Phi_Euler = Phi_Euler+1;
    end
end

Ramanujan_small = zeros(q,Phi_Euler);
for Number_Phi = 1:Phi_Euler
    Ramanujan_small(:,Number_Phi) = circshift(ramanujan,Number_Phi-1,2);
end
Number_Rep = P/q;
if q==1;
    Ramanujan_Space = repmat(Ramanujan_small,Number_Rep,1);
else
    Ramanujan_1 = repmat(Ramanujan_small,Number_Rep,1);

Ramanujan_Space = [Ramanujan_Space,Ramanujan_1];
end
end
for i=1:P
    for j=1:P
        if abs(Ramanujan_Space(i,j))<0.000000000001
            Ramanujan_Space(i,j) =0;
        end
    end
end

      





        
           
    
    
