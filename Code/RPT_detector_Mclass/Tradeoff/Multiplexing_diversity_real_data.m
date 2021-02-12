% Multiplexing Diversity concept and trade of between P_e and number of classes
%Probability of error vs. M
%Last edit March 2018
clear all
fs = 256;

Number_of_trials = 15;
Number_train = 3;
Harmonics = 3;
Number_test_class = Number_of_trials-Number_train;

% USCD dataset
Target_Freq = [9.25, 11.25, 9.75, 11.75, 10.25, 12.25,14.25, 10.75,12.75]; % 9 classes 
%Target_Set = [1  ,     2,  3  ,  4   ,   5  ,     6,    7,   8  ,     9];
T_set_length = length(Target_Freq);
Latency = fs/4;
L = 128;
Channels_set=[1,2,4,8];
SNR = 0.48;
for  N_ch= 1:4;
     Number_of_Channels = Channels_set(N_ch);
    counter_M = 1;
    [log_P_e_RPT] = RPT_Mutliplexing_diversity_real_data(Target_Freq,T_set_length,L,Number_test_class,Number_of_Channels,counter_M,fs,Number_of_trials,Number_train,Latency);    
    [log_P_e_CCA,log_P_e_IT_CCA] = CCA_Multiplexing_diversity_real_data(Target_Freq,T_set_length,L,Number_test_class,Number_of_Channels,counter_M,fs,Number_of_trials,Number_train,Harmonics,Latency);   
    figure(Number_of_Channels)
    plot(log2(2:T_set_length),-log_P_e_RPT/(L*SNR),'--s','Linewidth',2)
    hold on
    plot(log2(2:T_set_length),-log_P_e_CCA/(L*SNR),'--*','Linewidth',2)
    plot(log2(2:T_set_length),-log_P_e_IT_CCA/(L*SNR),'--d','Linewidth',2)
    grid on
    xlabel('log_2 M')
    ylabel('-log(P_e)/(L.SNR)')
    set_figure_size(900,700)
    legend('RPT','Standard CCA','IT CCA')
    figure(5)
    plot(log2(2:T_set_length),-log_P_e_RPT/(L*SNR),'--s','Linewidth',2)
    hold on
    grid on
    xlabel('log_2 M')
    ylabel('-log(P_e)/(L.SNR)')
    set_figure_size(900,700)
end
legend('N_c = 1','N_c = 2','N_c = 4','N_c = 8')



