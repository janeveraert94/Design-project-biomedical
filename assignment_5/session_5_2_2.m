% session_5_2_2
% November 2017

clear; close all; 
%clc; 

load epilepsy

%% Standardize if whished
standardize = 0;
if standardize
    data_train = zscore(data_train,[ ],2);
    data_test = zscore(data_test,[ ],2);
end

%% 1) Plot the training data in two adjacent subplots, one for each class
train_seiz_indx = find(labels_train == 1);
train_non_seiz_indx = find(labels_train == 0);

train_seiz_samples = data_train(train_seiz_indx,:);
train_non_seiz_samples = data_train(train_non_seiz_indx,:);

[linked_samples_train_seiz] = link_samples(train_seiz_samples);
[linked_samples_train_non_seiz] = link_samples(train_non_seiz_samples);

figure()
ax1 = subplot(211);
plot(linked_samples_train_seiz)
title('total of seizeure samples')
xlabel('samples')
ylabel('EEG [µV]')
axis tight
ax2 = subplot(212);
plot(linked_samples_train_non_seiz)
title('total of non seizeure samples')
xlabel('samples')
ylabel('EEG [µV]')
axis tight
linkaxes([ax1,ax2],'xy')

% Normalize or standardize????

%% 2) Calculate the signal energy for training and test segments.
for i = 1:size(data_train,1)
    sample = data_train(i,:);
    energy_sample = sum(sample.^2);
    energy_training_set(i,1) = energy_sample;
end

figure()
stem(train_seiz_indx,energy_training_set(train_seiz_indx)); hold on
stem(train_non_seiz_indx,energy_training_set(train_non_seiz_indx));hold off
axis tight
legend ('seizure samples','non seizure samples')
xlabel('sample number')
ylabel('signal energy [(µV)^2 * s ]')


%% 3) Calculate the frequency band energies based on the discrete wavelet decomposition
 %    Decompose each segment (training and test) into wavelet coefficients
 %    Calculate the energy contained in the coefficients for each band (A1, D1, D2,..., D5)
level_of_decomposition = 5;
[app_coef_train,det_coeff_train] = wavelet_dec(data_train,level_of_decomposition);
[app_coef_test,det_coeff_test] = wavelet_dec(data_test,level_of_decomposition);


%% 4) Display all the features as boxplots
 %    full energy of the signals and subband energies (7 subplots in total)

% calculate the energy of the subbands energies training set

% app_coef_train_energy
for i = 1:size(app_coef_train,1)
    sample = app_coef_train(i,:);
    energy_sample = sum(sample.^2);
    app_coef_training_energy(i,1) = energy_sample;
end

% det_coeff_train_energy
for i = 1:size(det_coeff_train,1)
    sample = det_coeff_train(i,:);
    for k = 1:size(sample,2)
        sample_coeff = cell2mat(sample(k));
        energy_sample_coeff = sum(sample_coeff.^2);
        det_coef_training_energy(i,k) = energy_sample_coeff;
    end
end
    
% plot
figure()
subplot(4,2,1)
boxplot(energy_training_set,labels_train);
title('original signal')
subplot(4,2,2)
boxplot(app_coef_training_energy,labels_train);
title('A5')
for k = 1:level_of_decomposition
    subplot(4,2,k+2);
    boxplot(det_coef_training_energy(:,k),labels_train);
    title(['D', num2str(level_of_decomposition-k+1)])
end





%% 5) Predict the test labels with the classify function
 %    using only the full signal energies
 %    using only the subband energies
 %    Note: would it make sense to combine both? Why (not)?

% calculation of the test samples signal energy
for i = 1:size(data_test,1)
    sample = data_test(i,:);
    energy_sample = sum(sample.^2);
    energy_test_set(i,1) = energy_sample;
end

% app_coef_test_energy
for i = 1:size(app_coef_test,1)
    sample = app_coef_test(i,:);
    energy_sample = sum(sample.^2);
    app_coef_test_energy(i,1) = energy_sample;
end
% det_coeff_test_energy
for i = 1:size(det_coeff_test,1)
    sample = det_coeff_test(i,:);
    for k = 1:size(sample,2)
        sample_coeff = cell2mat(sample(k));
        energy_sample_coeff = sum(sample_coeff.^2);
        det_coef_test_energy(i,k) = energy_sample_coeff;
    end
end

% classification LDA
LDA_signal_energy_label = classify(energy_test_set,energy_training_set,labels_train);
LDA_A5_energy_label = classify(app_coef_test_energy,app_coef_training_energy,labels_train);
for k = 1:size(det_coef_test_energy,2)
        LDA_Dcoef_energy = classify(det_coef_test_energy(:,k),det_coef_training_energy(:,k),labels_train);
        LDA_Dcoef_energy_label_all(:,k) = LDA_Dcoef_energy;
end

LDA_D5_D1_energy_label = classify([app_coef_test_energy det_coef_test_energy(:,1) det_coef_test_energy(:,2) det_coef_test_energy(:,3) det_coef_test_energy(:,4) det_coef_test_energy(:,5)],...
                [app_coef_training_energy det_coef_training_energy(:,1) det_coef_training_energy(:,2) det_coef_training_energy(:,3) det_coef_training_energy(:,4) det_coef_training_energy(:,5)],labels_train);

% LDA_all_energy_label = classify([energy_test_set app_coef_test_energy det_coef_test_energy(:,1) det_coef_test_energy(:,2) det_coef_test_energy(:,3) det_coef_test_energy(:,4) det_coef_test_energy(:,5)],...
%                 [energy_training_set app_coef_training_energy det_coef_training_energy(:,1) det_coef_training_energy(:,2) det_coef_training_energy(:,3) det_coef_training_energy(:,4) det_coef_training_energy(:,5)],labels_train);
            

% evaluate
check_LDA_signal_energy = LDA_signal_energy_label - labels_test;
number_of_mislabeled_LDA_signal_energy = nnz(check_LDA_signal_energy);
percentage_mislabeled_LDA_signal_energy = number_of_mislabeled_LDA_signal_energy/length(labels_test)*100;

check_LDA_A5_energy_label = LDA_A5_energy_label - labels_test;
number_of_mislabeled_LDA_A5_energy_label = nnz(check_LDA_A5_energy_label);
percentage_mislabeled_LDA_A5_energy_label = number_of_mislabeled_LDA_A5_energy_label/length(labels_test)*100;

for k = 1:size(LDA_Dcoef_energy_label_all,2)
    check_LDA_D_energy_label = LDA_Dcoef_energy_label_all(:,k) - labels_test;
    number_of_mislabeled_LDA_D_energy_label = nnz(check_LDA_D_energy_label);
    percentage_mislabeled_LDA_D_energy_label = number_of_mislabeled_LDA_D_energy_label/length(labels_test)*100;
    number_of_mislabeled_LDA_D_all_energy_label(1,k) = number_of_mislabeled_LDA_D_energy_label;
    percentage_mislabeled_LDA_D_all_energy_label(1,k) = percentage_mislabeled_LDA_D_energy_label;
end

check_LDA_D5_D1_energy_label = LDA_D5_D1_energy_label - labels_test;
number_of_mislabeled_LDA_D5_D1_energy_label = nnz(check_LDA_D5_D1_energy_label);
percentage_mislabeled_LDA_D5_D1_energy_label = number_of_mislabeled_LDA_D5_D1_energy_label/length(labels_test)*100;

% check_LDA_all_energy_label = LDA_all_energy_label - labels_test;
% number_of_mislabeled_LDA_all_energy_label = nnz(check_LDA_all_energy_label);
% percentage_mislabeled_LDA_all_energy_label = number_of_mislabeled_LDA_all_energy_label/length(labels_test)*100;


% plots
figure()
subplot(size(LDA_Dcoef_energy_label_all,2)+1,1,1)
plot_ft_with_false_labeling(app_coef_test_energy,LDA_A5_energy_label,...
                    labels_test);
ylabel(['A', num2str(level_of_decomposition)])
title('classification based on different coeff energy')
    legend ('labeled as seizure','labeled as non seizure','false labeled')


for k = 1:size(LDA_Dcoef_energy_label_all,2)
subplot(size(LDA_Dcoef_energy_label_all,2)+1,1,1+k)
plot_ft_with_false_labeling(det_coef_test_energy(:,k),LDA_Dcoef_energy_label_all(:,k),...
                    labels_test);
ylabel(['D', num2str(level_of_decomposition-k+1)])
end

figure()
plot_ft_with_false_labeling(energy_test_set,LDA_signal_energy_label,...
                    labels_test);
title('classification based on the signal energy')
ylabel('sample')
%% 6) Evaluate your classifiers (accuracy, sensitivity, specificity)
k = size(det_coef_test_energy,2);
sensitivity = zeros(k+3,1);
specificity = zeros(k+3,1);
accuracy = zeros(k+3,1);

for k = 1:size(det_coef_test_energy,2)
        [confusion_matrices{k+2},sensitivity(k+2),specificity(k+2),accuracy(k+2)] = confusion_matrix(LDA_Dcoef_energy_label_all(:,k),labels_test);
end
[confusion_matrices{2}, sensitivity(2), specificity(2), accuracy(2)] = confusion_matrix(LDA_A5_energy_label,labels_test);
[confusion_matrices{1}, sensitivity(1), specificity(1), accuracy(1)] = confusion_matrix(LDA_signal_energy_label,labels_test);
[confusion_matrices{k+3}, sensitivity(k+3), specificity(k+3), accuracy(k+3)] = confusion_matrix(LDA_D5_D1_energy_label,labels_test);
% [confusion_matrices{k+4}, sensitivity(k+4), specificity(k+4), accuracy(k+4)] = confusion_matrix(LDA_all_energy_label,labels_test);



% summary
Features = {'total signal';'A5';'D5';'D4';'D3';'D2';'D1';'D1, D2, D3, D4, D5, A5'};
table(Features,accuracy, sensitivity, specificity)


%% functions

% put samples next to each other
function [linked_samples] = link_samples(array_of_samples)
nb_of_samples = size(array_of_samples,1);
linked_samples = [];
for k = 1:nb_of_samples
    linked_samples = [linked_samples array_of_samples(k,:)];    
end
end

% wavelet decomp for array of signals
function [app_coef_array,det_coeff_array] = wavelet_dec(array_of_sigals,lvl_of_decomp)
for k = 1:size(array_of_sigals,1)
    sample = array_of_sigals(k,:);
    [C,L] = wavedec(sample,lvl_of_decomp,'db4');
    detcoef_array = detcoef(C,L,fliplr(1:lvl_of_decomp));
    A = appcoef(C,L,'db4',lvl_of_decomp);
    app_coef_array(k,:) = A;
    det_coeff_array(k,:) = detcoef_array;
end
end

% plot featurevector with indication of false labeling
function[] = plot_ft_with_false_labeling(features,label_vector,...
                    exact_label_vector)
    false_label_check = label_vector - exact_label_vector;
    ind_false_labels = find(false_label_check);
    
    label_seiz_indx = find(label_vector == 1);
    label_non_seiz_indx = find(label_vector == 0);

    % figures
    stem(label_seiz_indx,features(label_seiz_indx),'filled'); hold on
    stem(label_non_seiz_indx,features(label_non_seiz_indx),'filled');
    stem(ind_false_labels,features(ind_false_labels),'LineStyle','none',...
        'MarkerEdgeColor','black');hold off
    axis tight
%     legend ('labeled as seizure','labeled as non seizure','false labeled')
%     ylabel('signal energy')


end

% confusion matrix
function [C, sensitivity, specificity, accuracy] = confusion_matrix(pred_labels_bin, true_labels_bin)
    C = zeros(2,2);
    TP = 0;
    TN = 0;
    FP = 0;
    FN = 0;
    for i = 1:length(pred_labels_bin)
        if (pred_labels_bin(i) == 1 && true_labels_bin(i) == 1)
            TP = TP +1;
        elseif (pred_labels_bin(i) == 0 && true_labels_bin(i) == 0)
            TN = TN+1;
        elseif (pred_labels_bin(i) == 1 && true_labels_bin(i) == 0)
            FP = FP+1;
        elseif (pred_labels_bin(i) == 0 && true_labels_bin(i) == 1)
            FN = FN+1;
        end
    end
    sensitivity = TP/(TP+FN);
    specificity = TN/(TN + FP);
    accuracy = (TP+TN)/length(pred_labels_bin);
    FNF = FN/(TP + FN);
    FPF = FP/(TN + FP);
    C = [TP FP; FN TN];
end

