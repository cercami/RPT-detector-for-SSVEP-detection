% Multiplexing Diversity concept and trade of between P_e and number of classes
%Probability of error vs. M
%Last edit October 2018
clear all
Number_of_trials = 15;
Number_train = 3;
fs = 256;
Number_of_trials = 15;
Harmonics = 3;
Number_test_class = Number_of_trials-Number_train;
N_fold = 15;
M = 12;
Number_train = Number_of_trials - M;
Number_of_Classes = 9;
Num_sub = 10;
% USCD dataset
Target_Freq = [9.25, 11.25, 9.75, 11.75, 10.25, 12.25,14.25, 10.75,12.75]; % 9 classes 
%Target_Set = [1  ,     2,  3  ,  4   ,   5  ,     6,    7,   8  ,     9];
T_set_length = length(Target_Freq);
L = 256;
Channels_set=[1:8];
filenames_list = dir('directory');
n_files  = length(filenames_list);   
f_i = Target_Freq(1:Number_of_Classes);
Number_test =  Number_test_class * Number_of_Classes;  
SNR = 0.48;
        N_Sample = fs;
        Latency = fs/4;

for fold = 1:N_fold
    Indices = crossvalind('LeaveMout', Number_of_trials,M);
    for Sub_index = 1:n_files
        Number_of_Channels = 8;
        [Observation_Mat_i] = Load_USCD_Subject_Full(Sub_index,f_i,Number_of_Channels,Number_of_Classes);
        T=L/fs; 
        var_estimate_length = 38;
        L_tot = L + var_estimate_length + Latency;
        Observation_Mat = zeros(Number_of_Channels,L_tot,Number_of_trials,Number_of_Classes);

        for Class_Num = 1:Number_of_Classes
            for trials = 1:Number_of_trials;
                Observation_Mat_trial = zeros(Number_of_Channels,L_tot);
                Observation_Mat_trial(:,:) = Observation_Mat_i(Class_Num,1:Number_of_Channels,1:L_tot,trials);
                Observation_Mat(:,:,trials,Class_Num) = eegfilt(Observation_Mat_trial,fs,4,30,0,floor((L_tot/3)-1),0,'fir1');
                
            end
        end

        for  N_ch= 1:8;
             Number_of_Channels = Channels_set(N_ch);
            [P_e_sub_RPT] = RPT_Mutliplexing_diversity_kfold(Indices,Number_train,Target_Freq,T_set_length,L,Number_test_class,Number_of_Channels,fs,Observation_Mat,var_estimate_length);


            P_e_RPT(N_ch,:,Sub_index) = P_e_sub_RPT;

        end
    end
    

    hold_P_e_RPT(:,:,fold) = mean(P_e_RPT,3);

    % For statistical analysis 
    P_e_RPT_all_subjects(:,:,:,fold) = P_e_RPT;

    
end

log_P_e_RPT_avg = log(mean(hold_P_e_RPT,3));


%% For statistical analysis 
P_e_RPT_avg_all_subjects = mean(P_e_RPT_all_subjects,4);


P_e_RPT_avg = mean(P_e_RPT_avg_all_subjects,3);


P_e_RPT_all_elec = zeros(Number_of_Classes-1,Num_sub);
P_e_RPT_all_elec(:) = P_e_RPT_avg_all_subjects(4,:,:);
P_e_RPT_all_elec = P_e_RPT_all_elec';

%% Figures
for N_p = 1:8    

grid on
xlabel('log_2 M')
ylabel('-log(P_e)/(L.SNR)')
set_figure_size(900,700)
legend('RPT','Standard CCA','IT CCA','FBCCA (optimized)', 'FBCCA')
    %}
figure(5)
plot(log2(2:T_set_length),-log_P_e_RPT_avg(N_p,:)/(L*SNR),'--s','Linewidth',2)
hold on
grid on
xlabel('log_2 M')
ylabel('-log(P_e)/(L.SNR)')
set_figure_size(900,700)

legend('N_c = 1','N_c = 2','N_c = 4','N_c = 8')

end

