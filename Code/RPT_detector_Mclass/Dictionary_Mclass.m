% Building the Dictionary space fast.
% With considering the knowledge of supports in multiple possible
% hypotheses, this code generates the corresponding sub-space.
%% Dictionary Based method for GLRT
%% In this code, we are writing the Ramanujan basis in terms of Dictionary method
function [H,Matrix_D] = Dictionary_Mclass(L,Harmonics,Number_of_Classes)

%% To add another class, you have to add harmonics to the below array.


for h = 1:Number_of_Classes
    index = 1;
    N_Supp = length(Harmonics{h});
    Dictionary = [];
    Euler_Mat = [];
    Supprt = Harmonics{h};
    for u = 1:N_Supp
        p = Supprt(1,u);
        [Ramanujan_Space,~] = Ramanujan(p);
        [Euler] = Euler_fn(p);
        Branch = Ramanujan_Space(:,end-Euler+1:end);
        Branch_Extended = repmat(Branch,ceil(L/p),1);
        Sub_dictionary = Branch_Extended(1:L,:);
        Dictionary = [Dictionary,Sub_dictionary];
        Euler_Mat = [Euler_Mat,Euler];
        P_i = Supprt(1,u);
        rept = Euler;
        Vector_D(index:index+rept-1) = P_i^2;
        index = index + rept;
        end
    H{h} = Dictionary;
    Matrix_D{h}=diag(Vector_D);
    clear Dictionary
    clear Vector_D
end

        
        
    

        
        
    
    
    
    
    
    
    
    
    
    
   
    
    
    
    
   