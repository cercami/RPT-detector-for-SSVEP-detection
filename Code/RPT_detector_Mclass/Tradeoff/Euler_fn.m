function [Euler] = Euler_fn(divisors) 
n = size(divisors);
n = n(1,2);
for i=1:n
P = divisors(i);
Phi_Euler = 1;
    for coprime=2:P
        f_1 = factor(coprime);
        f_2 = factor(P);
        check = ismember(f_1,f_2); 
        check = sum(check);
        if check==0
           Phi_Euler = Phi_Euler+1;
        end
    end
    Euler(i) = Phi_Euler;
end

    