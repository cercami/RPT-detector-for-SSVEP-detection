%% Journal of Neural Engineering
% An extenstion to the Asilomar paper
% We address multiple classes, multiple channels
%% Define variable;
%% Last edit Feb 2019
function [Accuracy_RPT_fold] = RPT_kfold(Indices,Observation_Mat,Number_test_class,L,Latency,N_folds,var_estimate_length)

Number_of_Classes = 9;
Number_of_Channels = 8;
Number_of_Observation = Number_of_Classes*150;
Number_of_trials = 15;

fs=256; % Samping Frequency is 256 Hz.


%% USCD dataset
Target_Freq = [9.25, 11.25, 9.75, 11.75, 10.25, 12.25,14.25, 10.75, 12.75]; %9classes
Number_test =  Number_test_class * Number_of_Classes;
Number_train_class = Number_of_trials - Number_test_class; 

N_Sample = fs;

for Class_num = 1:Number_of_Classes
    T_i(Class_num) = round(fs/Target_Freq(Class_num));
    div_T{1,Class_num} = divsr(T_i(Class_num));
    f_i(Class_num) = Target_Freq(Class_num);
    Hp(:,Class_num)= Class_num*ones(Number_test_class,1);
end
True_Label = reshape(Hp,Number_test,1);
Counter_x = 1;

%% Cross validation
for fold = 1:N_folds
%% Divide data to Pre-stimulus and Post-Stimulus
    Test_set = (Indices==fold);
    Train_set = ~Test_set;

    Train_data = Observation_Mat(:,:,Train_set,:);
    Test_data = Observation_Mat(:,:,Test_set,:);

    Pre_stimulus_Train = Train_data(:,1:var_estimate_length,:,:);
    Post_stimulus_Train = Train_data(:,var_estimate_length+Latency+1:end,:,:);
    %Pre_stimulus_Test = Test_data(:,1:var_estimate_length,:,:);
    Post_stimulus_Test = Test_data(:,var_estimate_length+Latency+1:end,:,:); 

    %% Dictionary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P_max = L;
    [H,~] = Dictionary_Mclass(L,div_T,Number_of_Classes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for Class_num = 1:Number_of_Classes
        A{1,Class_num} = H{Class_num} * inv(H{Class_num}'*H{Class_num}) * H{Class_num}';
    end
    %% Covariance matrix estimation
  
        counter_test = 1;
        Concat_mean_Train =[];
        Concat_mean_Train_pre = [];
        Concat_mean_zero_Train_pre = [];

        for Class_num = 1:Number_of_Classes
            Mean_Train = mean(Post_stimulus_Train(:,:,:,Class_num),3);
            for Train_index = 1:Number_train_class
            Train_covariance = Post_stimulus_Train(:,:,Train_index,Class_num) - Mean_Train;
            Concat_mean_Train = [Concat_mean_Train,Train_covariance];
            Concat_mean_Train_pre = [Concat_mean_Train_pre,Pre_stimulus_Train(:,:,Train_index,Class_num)];
            % New method for covariance  % Added April 2018 %Uses pre stimulus
            Train_trial_mean = mean(Pre_stimulus_Train(:,:,Train_index,Class_num),2);
            Train_trial_mean_ext = repmat(Train_trial_mean,1,var_estimate_length);
            Train_mean_zero = Pre_stimulus_Train(:,:,Train_index,Class_num) - Train_trial_mean_ext;
            Concat_mean_zero_Train_pre = [Concat_mean_zero_Train_pre,Train_mean_zero];     
     
            end
        end
        
        Covariance_mat = cov(Concat_mean_Train_pre');
        Inv_Covarinace_mat = inv(Covariance_mat);
            
        for Class_ind = 1:Number_of_Classes
          
            for epoch = 1:Number_test_class;        


                y = zeros(Number_of_Channels,L);
                y(:,:) = Post_stimulus_Test(:,:,epoch,Class_ind);
                

                Y=y';
                tic
                for Class_num = 1:Number_of_Classes 
                    Suff_stat(Class_num) = 0.5*sum(dot((Y*Inv_Covarinace_mat)',Y'*A{1,Class_num}));
                end
                [~,Label(counter_test)] = max(Suff_stat);
                RPT_clock(counter_test) = toc;
                counter_test = counter_test + 1;
            end
        end
  %% Classification
    Correct_Label = (Label'==True_Label);
    TPN = sum(Correct_Label);
    Accuracy_RPT_fold(Counter_x) = (TPN)/(Number_test);
    RPT_clock_K(Counter_x) = mean(RPT_clock);
    Counter_x = Counter_x +1;
end % end to the fold loop
