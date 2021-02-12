function [P_e] = RPT_Mutliplexing_diversity_kfold(Indices,Number_train,Target_Freq,T_set_length,L,Number_test_class,Number_of_Channels,fs,Observation_Mat,var_estimate_length)
        counter_M = 1;
        for Number_of_Classes = 2:T_set_length
            f_i = Target_Freq(1:Number_of_Classes);
            Number_test =  Number_test_class * Number_of_Classes;       

            for Class_num = 1:Number_of_Classes
                T_i(Class_num) = round(fs/f_i(Class_num));
                div_T{1,Class_num} = divsr(T_i(Class_num));
                Hp(:,Class_num)= Class_num*ones(Number_test_class,1);
            end
            True_Label = reshape(Hp,Number_test,1);
            [H,~] = Dictionary_Mclass(L,div_T,Number_of_Classes);
            for Class_num = 1:Number_of_Classes
                 A{1,Class_num} = H{Class_num} * inv(H{Class_num}'*H{Class_num}) * H{Class_num}';
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
            N_Sample = fs;
            Latency = fs/4;

    
%% Divide data to Pre-stimulus and Post-Stimulus 
        Train_set = (Indices==1);
        Test_set = ~Train_set;
        Train_data = Observation_Mat(1:Number_of_Channels,:,Train_set,1:Number_of_Classes);
        Test_data = Observation_Mat(1:Number_of_Channels,:,Test_set,1:Number_of_Classes);
        Pre_stimulus_Train = Train_data(:,1:var_estimate_length,:,:);
        Post_stimulus_Train = Train_data(:,var_estimate_length+Latency+1:end,:,:);
        Post_stimulus_Test = Test_data(:,var_estimate_length+Latency+1:end,:,:);

%% Covariance Estimation
  
        counter_test = 1;
        Concat_mean_Train =[];
        Concat_mean_Train_pre = [];
        for Class_num = 1:Number_of_Classes
            Mean_Train = mean(Post_stimulus_Train(:,:,:,Class_num),3);
            for Train_index = 1:Number_train
            Train_covariance = Post_stimulus_Train(:,:,Train_index,Class_num) - Mean_Train;
            Concat_mean_Train = [Concat_mean_Train,Train_covariance];
            Concat_mean_Train_pre = [Concat_mean_Train_pre,Pre_stimulus_Train(:,:,Train_index,Class_num)];
            end
        end
        Covariance_mat = cov(Concat_mean_Train_pre');
        Inv_Covarinace_mat = inv(Covariance_mat);
 
%% Classification
        for Class_ind = 1:Number_of_Classes          
            for epoch = 1:Number_test_class; 
                y = zeros(Number_of_Channels,L);
                y(:,:) = Post_stimulus_Test(:,:,epoch,Class_ind);

                Y=y';
                for Class_num = 1:Number_of_Classes 
                    Suff_stat(Class_num) = 0.5*sum(dot((Y*Inv_Covarinace_mat)',Y'*A{1,Class_num}));
                end
                [~,Label(counter_test)] = max(Suff_stat);
                counter_test = counter_test + 1;
            end
        end
        Correct_Label = (Label'==True_Label);
        TPN = sum(Correct_Label);
        Errors = (Number_test - TPN)/Number_test_class;
        P_e(counter_M) = (1/Number_of_Classes)*Errors;
        counter_M = counter_M + 1;
        clear Label
    end
        
   
    

