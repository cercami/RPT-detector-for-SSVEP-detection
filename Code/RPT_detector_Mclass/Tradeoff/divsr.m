function [divisor] = divsr (P)
counter = 1;
for d=1:P
    if mod(P,d)==0
        divisor(1,counter)=d;
        counter=counter+1;
    end
end

