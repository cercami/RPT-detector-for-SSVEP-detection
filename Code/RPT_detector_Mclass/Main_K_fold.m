%% Main
% This code is for M-fold validation of RPT, CCA IT and CCA
clear all
Number_of_trials = 15;
N_folds = 5;
Number_test_class = Number_of_trials/N_folds;
Number_train_class = Number_of_trials - Number_test_class;
num_sub = 10;
Number_of_Classes = 9;
Number_of_Channels = 8;
fs = 256;
Indices = crossvalind('kfold', Number_of_trials,N_folds);
Latency = fs/4;
% Target frequencies
Target_Freq = [9.25, 11.25, 9.75, 11.75, 10.25, 12.25,14.25, 10.75, 12.75]; %9classes

Number_test =  Number_test_class * Number_of_Classes;

%% Load Data
filenames_list = dir('input_directory');
n_files  = length(filenames_list);
N_Sample = fs;

for Sub_index = 1:n_files
    
    %% Load data
    [Observation_Mat_i] = Load_USCD_Subject_Full(Sub_index,Target_Freq,Number_of_Channels,Number_of_Classes);
    
    counter_t = 1;
    for K = 0.5:0.5:3.5;
        L = K*N_Sample;
        T=L/fs; 
        var_estimate_length = 38;
        L_tot = L + var_estimate_length+Latency;
        Observation_Mat = zeros(Number_of_Channels,L_tot,Number_of_trials,Number_of_Classes);
        %% Preprocessing
        for Class_Num = 1:Number_of_Classes
            for trials = 1:Number_of_trials;
                Observation_Mat_trial = zeros(Number_of_Channels,L_tot);
                Observation_Mat_trial(:,:) = Observation_Mat_i(Class_Num,1:Number_of_Channels,1:L_tot,trials);
                Observation_Mat(:,:,trials,Class_Num) = eegfilt(Observation_Mat_trial,fs,4,30,0,floor((L_tot/3)-1),0,'fir1');
            end
        end
    
        %RPT    
        [Accuracy_RPT_fold] = RPT_kfold(Indices,Observation_Mat,Number_test_class,L,Latency,N_folds,var_estimate_length);
        Accuracy_RPT_mat(Sub_index,counter_t,:) = Accuracy_RPT_fold;    

   
        %IT CCA
        [Accuracy_ITCCA_fold] = CCAIT_kfold(Indices,Observation_Mat,Number_test_class,L,T,Latency,N_folds,var_estimate_length);
        Accuracy_ITCCA_mat(Sub_index,counter_t,:) = Accuracy_ITCCA_fold;

        % Standard CCA
        [Accuracy_CCA_fold] = Standard_CCA_kfold(Indices,Observation_Mat,Number_test_class,L,T,Latency,N_folds,var_estimate_length);
        Accuracy_CCA_mat(Sub_index,counter_t,:) = Accuracy_CCA_fold;
        counter_t = counter_t+1;
        end
    % FBCCA optimized and non-optimized
    [Accuracy_FBCCAO_time,Accuracy_FBCCA_time] = FBCCA_optimized_kfold(Indices,Observation_Mat_i,Number_test_class,Latency,N_folds);
    Accuracy_FBCCA_optimized_mat(Sub_index,:,:) = Accuracy_FBCCAO_time;
    Accuracy_FBCCA_mat(Sub_index,:,:) = Accuracy_FBCCA_time;

    % PSDA
    %[Accuracy_PSDA_time] = PSDA_kfold(Indices,Observation_Mat_i,Number_test_class,Latency,N_folds);
    %Accuracy_PSDA_mat(Sub_index,:,:) = Accuracy_PSDA_time;
        
end

K = 0.5:0.5:3.5;
K_mat = repmat(K,num_sub,1);
%RPT
Accuracy_avg_sub_RPT = mean(Accuracy_RPT_mat,3);
Accuracy_avg_RPT = mean(Accuracy_avg_sub_RPT);

ITR_log_RPT = Accuracy_avg_sub_RPT.*log2(Accuracy_avg_sub_RPT) + (1-Accuracy_avg_sub_RPT).*log2((1-Accuracy_avg_sub_RPT)/(Number_of_Classes-1));
ITR_log_RPT(isnan(ITR_log_RPT))=0;
ITR_Final_RPT_sub = (60./K_mat).*(log2(Number_of_Classes) + ITR_log_RPT);

% IT CCA
Accuracy_avg_sub_ITCCA = mean(Accuracy_ITCCA_mat,3);
Accuracy_avg_ITCCA = mean(Accuracy_avg_sub_ITCCA);
ITR_log_ITCCA = Accuracy_avg_sub_ITCCA.*log2(Accuracy_avg_sub_ITCCA) + (1-Accuracy_avg_sub_ITCCA).*log2((1-Accuracy_avg_sub_ITCCA)/(Number_of_Classes-1));
ITR_log_ITCCA(isnan(ITR_log_ITCCA))=0;
ITR_Final_ITCCA_sub = (60./K_mat).*(log2(Number_of_Classes) + ITR_log_ITCCA);

% Standard CCA
Accuracy_avg_sub_CCA = mean(Accuracy_CCA_mat,3);
Accuracy_avg_CCA = mean(Accuracy_avg_sub_CCA);
ITR_log_CCA = Accuracy_avg_sub_CCA.*log2(Accuracy_avg_sub_CCA) + (1-Accuracy_avg_sub_CCA).*log2((1-Accuracy_avg_sub_CCA)/(Number_of_Classes-1));
ITR_log_CCA(isnan(ITR_log_CCA))=0;
ITR_Final_CCA_sub = (60./K_mat).*(log2(Number_of_Classes) + ITR_log_CCA);

% FBCCA Optimized
Accuracy_avg_sub_FBCCA_optimized = mean(Accuracy_FBCCA_optimized_mat,3);
Accuracy_avg_FBCCA_optimized = mean(Accuracy_avg_sub_FBCCA_optimized);
ITR_log_FBCCA_optimized = Accuracy_avg_sub_FBCCA_optimized.*log2(Accuracy_avg_sub_FBCCA_optimized) + (1-Accuracy_avg_sub_FBCCA_optimized).*log2((1-Accuracy_avg_sub_FBCCA_optimized)/(Number_of_Classes-1));
ITR_log_FBCCA_optimized(isnan(ITR_log_FBCCA_optimized))=0;
ITR_Final_FBCCA_sub_optimized = (60./K_mat).*(log2(Number_of_Classes) + ITR_log_FBCCA_optimized);

% FBCCA
Accuracy_avg_sub_FBCCA = mean(Accuracy_FBCCA_mat,3);
Accuracy_avg_FBCCA = mean(Accuracy_avg_sub_FBCCA);
ITR_log_FBCCA = Accuracy_avg_sub_FBCCA.*log2(Accuracy_avg_sub_FBCCA) + (1-Accuracy_avg_sub_FBCCA).*log2((1-Accuracy_avg_sub_FBCCA)/(Number_of_Classes-1));
ITR_log_FBCCA(isnan(ITR_log_FBCCA))=0;
ITR_Final_FBCCA_sub = (60./K_mat).*(log2(Number_of_Classes) + ITR_log_FBCCA);

%PSDA
%K = 1:1:3;
%K_mat = repmat(K,num_sub,1);

%Accuracy_avg_sub_PSDA = mean(Accuracy_PSDA_mat,3);
%Accuracy_avg_PSDA = mean(Accuracy_avg_sub_PSDA);
%ITR_log_PSDA = Accuracy_avg_sub_PSDA.*log2(Accuracy_avg_sub_PSDA) + (1-Accuracy_avg_sub_PSDA).*log2((1-Accuracy_avg_sub_PSDA)/(Number_of_Classes-1));
%ITR_log_PSDA(isnan(ITR_log_PSDA))=0;
%ITR_Final_PSDA_sub = (60./K_mat).*(log2(Number_of_Classes) + ITR_log_PSDA);

%% Figures 
K = 0.5:0.5:3.5;
figure(1)
plot(K,Accuracy_avg_RPT,'--s','Color',[0 0.45 0.74],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[0 0.45 0.74],'MarkerFaceColor',[0 0.45 0.74])
hold on
plot(K,Accuracy_avg_ITCCA,':d','Color',[1 0 0],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1,0,0])
plot(K,Accuracy_avg_CCA,'-.o','Color',[0.47 0.67 0.19],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[0.47 0.67 0.19],'MarkerFaceColor',[0.47 0.67 0.19])
plot(K,Accuracy_avg_FBCCA_optimized,':*','Color',[0.49 0.18 0.56],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[0.49 0.18 0.56],'MarkerFaceColor',[0.49 0.18 0.56])
plot(K,Accuracy_avg_FBCCA,'--v','Color',[0 0 0],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0,0,0])
%K = 1:1:3;
%plot(K,Accuracy_avg_PSDA,'--^','Color',[0.93 0.69 0.13],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[0.93 0.69 0.13],'MarkerFaceColor',[0.93 0.69 0.13])

xlabel('Data length (Seconds)')
ylabel('Accuracy')
legend('RPT detector','IT CCA','Standard CCA','FBCCA (optimized)','FBCCA')
grid on
set_figure_size(900,700)

K = 0.5:0.5:3.5;
figure(2)
errorbar(K,mean(Accuracy_avg_sub_RPT),sqrt(var(Accuracy_avg_sub_RPT))/sqrt(num_sub),'--s','Color',[0 0.45 0.74],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[0 0.45 0.74],'MarkerFaceColor',[0 0.45 0.74])
hold on
errorbar(K,mean(Accuracy_avg_sub_ITCCA),sqrt(var(Accuracy_avg_sub_ITCCA))/sqrt(num_sub),':d','Color',[1 0 0],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1,0,0])
errorbar(K,mean(Accuracy_avg_sub_CCA),sqrt(var(Accuracy_avg_sub_CCA))/sqrt(num_sub),'-.o','Color',[0.47 0.67 0.19],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[0.47 0.67 0.19],'MarkerFaceColor',[0.47 0.67 0.19])
errorbar(K,mean(Accuracy_avg_sub_FBCCA_optimized),sqrt(var(Accuracy_avg_sub_FBCCA_optimized))/sqrt(num_sub),':*','Color',[0.49 0.18 0.56],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[0.49 0.18 0.56],'MarkerFaceColor',[0.49 0.18 0.56])
errorbar(K,mean(Accuracy_avg_sub_FBCCA),sqrt(var(Accuracy_avg_sub_FBCCA))/sqrt(num_sub),'--v','Color',[0 0 0],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0,0,0])
%K = 1:1:3;
%errorbar(K,mean(Accuracy_avg_sub_PSDA),sqrt(var(Accuracy_avg_sub_PSDA))/sqrt(num_sub),'--^','Color',[0.93 0.69 0.13],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[0.93 0.69 0.13],'MarkerFaceColor',[0.93 0.69 0.13])
xlabel('Data length (Seconds)')
ylabel('Accuracy')
legend('RPT detector','IT CCA','Standard CCA','FBCCA (optimized)','FBCCA')
grid on
set_figure_size(900,700)

K = 0.5:0.5:3.5;
figure(3)
errorbar(K,mean(ITR_Final_RPT_sub),sqrt(var(ITR_Final_RPT_sub))/sqrt(num_sub),'--s','Color',[0 0.45 0.74],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[0 0.45 0.74],'MarkerFaceColor',[0 0.45 0.74])
hold on
errorbar(K,mean(ITR_Final_ITCCA_sub),sqrt(var(ITR_Final_ITCCA_sub))/sqrt(num_sub),':d','Color',[1 0 0],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1,0,0])
errorbar(K,mean(ITR_Final_CCA_sub),sqrt(var(ITR_Final_CCA_sub))/sqrt(num_sub),'-.o','Color',[0.47 0.67 0.19],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[0.47 0.67 0.19],'MarkerFaceColor',[0.47 0.67 0.19])
errorbar(K,mean(ITR_Final_FBCCA_sub_optimized),sqrt(var(ITR_Final_FBCCA_sub_optimized))/sqrt(num_sub),':*','Color',[0.49 0.18 0.56],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[0.49 0.18 0.56],'MarkerFaceColor',[0.49 0.18 0.56])
errorbar(K,mean(ITR_Final_FBCCA_sub),sqrt(var(ITR_Final_FBCCA_sub))/sqrt(num_sub),'--v','Color',[0 0 0],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0,0,0])
K = 1:1:3;
%errorbar(K,mean(ITR_Final_PSDA_sub),sqrt(var(ITR_Final_PSDA_sub))/sqrt(num_sub),'--^','Color',[0.93 0.69 0.13],'linewidth',2,'MarkerSize',8,'MarkerEdgeColor',[0.93 0.69 0.13],'MarkerFaceColor',[0.93 0.69 0.13])

xlabel('Data length (Seconds)')
ylabel('ITR (bits/min)')
legend('RPT detector','IT CCA','Standard CCA','FBCCA (optimized)','FBCCA')
grid on
set_figure_size(900,700)

K = 0.5:0.5:3.5;
figure(4)
plot(K,mean(ITR_Final_RPT_sub),'--s','Color',[0 0.45 0.74],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[0 0.45 0.74],'MarkerFaceColor',[0 0.45 0.74])
hold on
plot(K,mean(ITR_Final_ITCCA_sub),':d','Color',[1 0 0],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1,0,0])
plot(K,mean(ITR_Final_CCA_sub),'-.o','Color',[0.47 0.67 0.19],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[0.47 0.67 0.19],'MarkerFaceColor',[0.47 0.67 0.19])
plot(K,mean(ITR_Final_FBCCA_sub_optimized),':*','Color',[0.49 0.18 0.56],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[0.49 0.18 0.56],'MarkerFaceColor',[0.49 0.18 0.56])
plot(K,mean(ITR_Final_FBCCA_sub),'--v','Color',[0 0 0],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0,0,0])

% K= 1:1:3;
%plot(K,mean(ITR_Final_PSDA_sub),'--^','Color',[0.93 0.69 0.13],'linewidth',3,'MarkerSize',8,'MarkerEdgeColor',[0.93 0.69 0.13],'MarkerFaceColor',[0.93 0.69 0.13])
xlabel('Data length (Seconds)')
ylabel('ITR (bits/min)')
legend('RPT detector','IT CCA','Standard CCA','FBCCA (optimized)','FBCCA')
grid on
set_figure_size(900,700)
[rhos_acc,rhos_ITR,h_acc,h_ITR] = Statistical_Analysis(Accuracy_avg_sub_RPT,ITR_Final_RPT_sub,Accuracy_avg_sub_ITCCA,ITR_Final_ITCCA_sub,Accuracy_avg_sub_CCA,ITR_Final_CCA_sub,Accuracy_avg_sub_FBCCA,ITR_Final_FBCCA_sub);
