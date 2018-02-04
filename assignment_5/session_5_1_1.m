% session_5_1_1.m
% November 2017
%%

clear; close all; clc; 
addpath ./FastICA_25/
rng('default');

%% Data
fs = 500;
t = (0:4999 )/500; 
freq = linspace(0,fs,length(t));

load icadata.mat
ecg = Sources(1,:);
reference = Sources(3,:); 

all_sources = Sources;
for k = 1:size(all_sources,1)
    sources_norm(k,:) = all_sources(k,:)-mean(all_sources(k,:));
    sources_norm(k,:) = all_sources(k,:)./std(all_sources(k,:));
end

%mix the sources
X = A*Sources;

%display the mixtures
figure()
for i = 1:size(X,1)
    subplot(size(X,1),1,i);
    plot(t, X(i,:));
    title(['mixture signal #', num2str(i)]);
    xlabel('Time (s)');
    ylabel('ECG [mV]');
end

%% Source separation

% 1) Apply fastICA on the mixtures in X to obtain the originals, save the components in a variable y 

[y] = fastica(X);
        
        
 %% Display the source estimates
 
%  figure()
% for i = 1:size(y,1)
%     subplot(size(y,1),1,i);
%     plot(t, y(i,:));
%     title(['icasig signal #', num2str(i)]);
%     xlabel('Time (s)');
% end
 

%% Normalization, pairing and performance

% 2) Normalize the y-estimates: set mean value to 0 and variance to 1. Do you know why?

for k = 1:size(y,1)
    y_norm(k,:) = y(k,:)-mean(y(k,:));
    y_norm(k,:) = y_norm(k,:)./std(y_norm(k,:));
end

% 3) Automatically match each Source with (one of) its best estimate(s))
%    Use squared or absolute correlation as matching criterion

[y1_resort_resigned,corresponding_sources] =...
    link_ICA_to_sources(y_norm,sources_norm);

figure('name','X1')
nb_of_ica = size(y1_resort_resigned,1);
for i = 1:nb_of_ica
   subplot(nb_of_ica,2,2*i-1)
   plot(t,sources_norm(corresponding_sources(i),:))
   title(strcat('source',num2str(corresponding_sources(i))))
end
xlabel('time [s]')
for i = 1:nb_of_ica
   subplot(nb_of_ica,2,2*i)
   plot(t,y1_resort_resigned(i,:))
   title(strcat('ICA signal',num2str(i)))
   
end
xlabel('time [s]')

% 4) calculate the Root Mean Squared Error (RMSE) between matched Source and estimate
%    y(t). Pay attention that the signals are zero mean, standardized and that they have the same sign. 
%    That is, use lower RMSE value (e.g. RMSE (Source_i, y_j ) and RMSE (Source_i, -y_j))
for i = 1:size(y1_resort_resigned,1)
    RMSE = sqrt(mean((sources_norm(i,:) - y1_resort_resigned(i,:)).^2));
    RMSE_all(i) = RMSE;
end
figure()
plot(sources_norm(4,:));hold on
plot(y1_resort_resigned(4,:)); hold off
ICAS = {'IC 1';'IC 2';'IC 3';'IC 4'};
RMSE = RMSE_all';
T = table(ICAS,RMSE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Number of components: other mixing matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X2 = A2*Sources;
X3 = A2*Sources + 0.1*rand(6,length(Sources));

figure()
for i = 1:size(X2,1)
    subplot(size(X2,1),1,i);
    plot(t, X2(i,:));
    title(['mixture signal  #', num2str(i)]);
    ylabel('ECG [mV]');
end
xlabel('Time (s)');
% 5) Apply fastICA for mixtures X2 and X3. Perform the normalization of the components and match them to the Sources.
%    What do you observe? How does it work or why doesn't it?

[y2] = fastica(X2);
[y3] = fastica(X3);

% for signal y2
for k = 1:size(y2,1)
    y2_norm(k,:) = y2(k,:)-mean(y2(k,:));
    y2_norm(k,:) = y2_norm(k,:)./std(y2_norm(k,:));
end

[ICA_2_resort_resigned,corresponding_sources] =...
    link_ICA_to_sources(y2_norm,sources_norm);

figure('name','X2')
nb_of_ica = size(ICA_2_resort_resigned,1);
for i = 1:nb_of_ica
   subplot(nb_of_ica,2,2*i-1)
   plot(t,sources_norm(corresponding_sources(i),:))
   title(strcat('source',num2str(corresponding_sources(i))))
end
xlabel('time [s]')
for i = 1:nb_of_ica
   subplot(nb_of_ica,2,2*i)
   plot(t,ICA_2_resort_resigned(i,:))
   title(strcat('ICA signal',num2str(i)))     
end
xlabel('time [s]')

% rmse for x2
for i = 1:size(ICA_2_resort_resigned,1)
    RMSE_2 = sqrt(mean((sources_norm(i,:) - ICA_2_resort_resigned(i,:)).^2));
    RMSE_2_all(i) = RMSE_2;
end
RMSE_2 = RMSE_2_all';
ICAS = {'IC 1';'IC 2';'IC 3';'IC 4'};
T = table(ICAS,RMSE_2)

% for signal y3
figure()
for i = 1:size(X3,1)
    subplot(size(X3,1),1,i);
    plot(t, X3(i,:));
    title(['mixture signal  #', num2str(i)]);
    ylabel('ECG [mV]');
end
xlabel('Time (s)');

for k = 1:size(y3,1)
    y3_norm(k,:) = y3(k,:)-mean(y3(k,:));
    y3_norm(k,:) = y3_norm(k,:)./std(y3_norm(k,:));
end

[ICA_3_resort_resigned,corresponding_sources] =...
    link_ICA_to_sources(y3_norm,sources_norm);

figure('name','X3')
nb_of_ica = size(ICA_3_resort_resigned,1);
for i = 1:nb_of_ica
   subplot(nb_of_ica,2,2*i-1)
   plot(t,sources_norm(corresponding_sources(i),:))
   title(strcat('source',num2str(corresponding_sources(i))))
end
xlabel('time [s]')
for i = 1:nb_of_ica
   subplot(nb_of_ica,2,2*i)
   plot(t,ICA_3_resort_resigned(i,:))
   title(strcat('ICA signal',num2str(i)))     
end
xlabel('time [s]')
% rmse for x2
for i = 1:size(ICA_3_resort_resigned,1)
    RMSE_3 = sqrt(mean((sources_norm(corresponding_sources(i),:)...
                - ICA_3_resort_resigned(i,:)).^2));
    RMSE_3_all(i) = RMSE_3;
end
RMSE_3 = RMSE_3_all';
ICAS = {'IC 1';'IC 2';'IC 3';'IC 4';'IC 5';'IC 6'};
T = table(ICAS,RMSE_3)

%% Artefact removal
% 6) Use fastICA to remove the "sawtooth" signal and reconstruct the mixtures signal. 
%    The sawtooth signal is saved in the variable reference
%    so USE this variable 'reference' to find the sawtooth component and
%    change the estimated mixing matrix A in order to remove it from the mixtures
%    hint: x=A*y (should y be normalized here, why (not)?)

clear Sources X2 X3 X4 %we don't use the original sources anymore 
% use the original mixture X and find the components with fastICA  

[y,~,W] = fastica(X);

for k = 1:size(y,1)
    y_norm(k,:) = y(k,:)-mean(y(k,:));
    y_norm(k,:) = y_norm(k,:)./std(y_norm(k,:));
end

% sort the ysignals and the W matrix corresponding to the sources
[y1_resort_resigned,corresponding_sources,W_resort] =...
    link_ICA_to_sources_and_W(y_norm,sources_norm,W);

figure('name','X1')
nb_of_ica = size(y1_resort_resigned,1);
for i = 1:nb_of_ica
   subplot(nb_of_ica,2,2*i-1)
   plot(t,sources_norm(corresponding_sources(i),:))
   title(strcat('source',num2str(corresponding_sources(i))))
end
xlabel('time [s]')
for i = 1:nb_of_ica
   subplot(nb_of_ica,2,2*i)
   plot(t,y1_resort_resigned(i,:))
   title(strcat('ICA signal',num2str(i)))     
end
xlabel('time [s]')


% make subplots of the original mixtures next to the mixtures from which the
% sawtooth has been removed

W_resort_inv = inv(W_resort);
W_resort_inv(:,3) = 0;
X_without_saw = W_resort_inv*y1_resort_resigned;
[y_without_saw,~,W_whitout_saw] = fastica(X_without_saw);

for k = 1:size(y_without_saw,1)
    y_without_saw_norm(k,:) = y_without_saw(k,:)-mean(y_without_saw(k,:));
    y_without_saw_norm(k,:) = y_without_saw_norm(k,:)./std(y_without_saw_norm(k,:));
end

[y_whitout_saw_resort_resigned,corresponding_sources,W_whitout_saw_resort] =...
    link_ICA_to_sources_and_W(y_without_saw_norm,sources_norm,W_whitout_saw);

figure('name','X1_without_saw')
nb_of_ica = size(y_whitout_saw_resort_resigned,1);
for i = 1:nb_of_ica
   subplot(nb_of_ica,2,2*i-1)
   plot(t,sources_norm(corresponding_sources(i),:))
   title(strcat('source',num2str(corresponding_sources(i))))
   subplot(nb_of_ica,2,2*i)
   plot(t,y_whitout_saw_resort_resigned(i,:))
   title(strcat('ICA signal',num2str(i)))     
end

figure()
subplot(211)
plot(t(500:1500),-1*X(3,500:1500)); hold on
plot(t(500:1500),X_without_saw(3,500:1500)); hold off
axis tight
xlabel('time [s]')
ylabel('ECG [µV]')
title('Mixed signal (time domain)')
legend('original','without sawtooth')
subplot(212)
plot(t(500:1500),-1*X(3,500:1500)-X_without_saw(3,500:1500));
axis tight
xlabel('time [s]')
ylabel('delta(ECG) [µV]')
title('Difference original-sawtooth removed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional filtering on signal X (sawtooth already removed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7) Choose one of the mixtures after the sawtooth has been filtered out.
%    Use appropriate filter techniques to clean it up and plot the result together with the original ecg signal (ecg).
%    Add plots to justify your filter choices (e.g. frequency spectrum, ...)

% fourier transform of the signal before filtering
for i = 1:size(X_without_saw,1)
    X_without_saw_ft = fft(X_without_saw(i,:));
    maxft = max(abs(X_without_saw_ft));
    X_without_sawspec = 20*log10(abs(X_without_saw_ft)/maxft);
    X_without_saw_spec(i,:) = X_without_sawspec;
end
% fourier transform of sources
for i = 1:size(sources_norm,1)
    sources_norm_ft = fft(sources_norm(i,:));
    maxft = max(abs(sources_norm_ft));
    sources_normspec = 20*log10(abs(sources_norm_ft)/maxft);
    sources_norm_spec(i,:) = sources_normspec;
end

% Plot the sources - time domain
figure()
for i = 1:size(sources_norm,1)
    subplot(size(sources_norm,1),2,2*i-1)
    plot(t,sources_norm(i,:))
    title(['Source #', num2str(i), '- time domain']);
end
xlabel('time [s]')
% Plot the sources - fr domain
for i = 1:size(sources_norm,1)
    subplot(size(sources_norm_spec,1),2,2*i)
    plot(freq(1:length(freq)/2),(sources_norm_spec(i,1:length(freq)/2)))
    axis tight
    title(['Source #', num2str(i), '- freq domain']);
    ylabel('Magn resp [dB]')

end
xlabel('Frequenty [Hz]')

% Plot the signals to be filtered - time domain
figure()
for i = 1:size(X_without_saw,1)
    subplot(size(X_without_saw,1),2,2*i-1)
    plot(t,X_without_saw(i,:))
end
xlabel('time [s]')
% Plot the signals to be filtered - fr domain
for i = 1:size(X_without_saw,1)
    subplot(size(X_without_saw,1),2,2*i)
    plot(freq(1:length(freq)/2),(X_without_saw_spec(i,1:length(freq)/2)))
    axis tight
end
ylabel('Magn resp [dB]')


%% Removal of the powerline interference -> Notch filter
% -> Notch voor die ene freq = 50Hz
f_0 = 50;
theta_0 = 2*pi*f_0/fs;
p = 0.99;
[b] = [0 0 1 -2*cos(theta_0) 1];
[a] = [1 -2*p*cos(theta_0) p^2 0 0];
h = fvtool(b,a);
h.Analysis = 'freq';
h.Fs = 250;
h.FrequencyRange='[0, Fs/2)';
grp_delay_notch = round(mean(grpdelay(b,a)))+1;
for i = 1:size(X_without_saw,1)
    output_notch = filter(b,a,X_without_saw(i,:));
    X_after_notch(i,:) = output_notch;
end

% freq domain after notch
for i = 1:size(X_after_notch,1)
    X_after_notch_ft = fft(X_after_notch(i,:));
    maxft = max(abs(X_after_notch_ft));
    X_after_notchspec = 20*log10(abs(X_after_notch_ft)/maxft);
    X_after_notch_spec(i,:) = X_after_notchspec;
end

% Plot the signals to be filtered - time domain
figure()
for i = 1:size(X_after_notch,1)
    subplot(size(X_after_notch,1),2,2*i-1)
    plot(t,X_after_notch(i,:))
end
% Plot the signals to be filtered - freq domain
for i = 1:size(X_after_notch_spec,1)
    subplot(size(X_after_notch,1),2,2*i)
    plot(freq(1:length(freq)/2),(X_after_notch_spec(i,1:length(freq)/2)))
    axis tight
end

figure()
subplot(221)
plot(t(500:950),X_without_saw(1,500:950))
axis tight
xlabel('time [s]')
title('Mixed signal (time domain)')
subplot(223)
plot(freq(2:length(freq)/2),(X_without_saw_spec(1,2:length(freq)/2)))
axis tight
xlabel('freq [Hz]')
ylabel('Magn resp [dB]')
title('Mixed signal (freq domain)')

subplot(222)
plot(t(500:950),X_after_notch(1,500:950))
axis tight
xlabel('time [s]')
title('Mixed signal after notch (time domain)')
subplot(224)
plot(freq(2:length(freq)/2),(X_after_notch_spec(1,2:length(freq)/2)))
axis tight
xlabel('freq [Hz]')
ylabel('Magn resp [dB]')
title('Mixed signal after notch (freq domain)')


%% -> LP voor die random HF noise

% Butterworth fc = 50 and n = 3
% dit heb ik gwn gekozen omdat dat in de tweede oefensessie het beste
% uitkwam -> nog checken hier!!!

% f_c = 50, n = 3
fc = 50;
co = fc/(fs/2);
[b3,a3] = butter(3,co);
for i = 1:size(X_after_notch,1)
    output_BW = filter(b3,a3,X_after_notch(i,:));
    X_after_all_filters(i,:) = output_BW;
end
% freq domain after BW
for i = 1:size(X_after_all_filters,1)
    X_after_all_filters_ft = fft(X_after_all_filters(i,:));
    maxft = max(abs(X_after_all_filters_ft));
    X_after_all_filtersspec = 20*log10(abs(X_after_all_filters_ft)/maxft);
    X_after_all_filters_spec(i,:) = X_after_all_filtersspec;
end

% Plot the signals to be filtered - time domain
figure()
for i = 1:size(X_after_all_filters,1)
    subplot(size(X_after_all_filters,1),2,2*i-1)
    plot(t,X_after_all_filters(i,:))
end
% Plot the signals to be filtered - freq domain
for i = 1:size(X_after_all_filters_spec,1)
    subplot(size(X_after_all_filters_spec,1),2,2*i)
    plot(freq(1:length(freq)/2),X_after_all_filters_spec(i,1:length(freq)/2))
    axis tight
end
grp_delay_butter = round(mean(grpdelay(b3,a3)))+1;

figure()
subplot(221)
plot(t(500:950),X_after_notch(1,500:950))
axis tight
xlabel('time [s]')
title('Mixed signal after notch (time domain)')
subplot(223)
plot(freq(1:length(freq)/2),(X_after_notch_spec(1,1:length(freq)/2)))
axis([0 250 -90 0])
xlabel('freq [Hz]')
ylabel('Magn resp [dB]')
title('Mixed signal after notch (freq domain)')

subplot(222)
plot(t(500:950),X_after_all_filters(i,500:950))
axis tight
xlabel('time [s]')
title('Mixed signal after notch & LP(time domain)')
subplot(224)
plot(freq(1:length(freq)/2),X_after_all_filters_spec(i,1:length(freq)/2))
axis([0 250 -90 0])
xlabel('freq [Hz]')
ylabel('Magn resp [dB]')
title('Mixed signal after notch & LP(freq domain)')

%% ICA on the mixture signals after all the filters

[y_after_all_filters,~,W_after_all_filters] = fastica(X_after_all_filters);

for k = 1:size(y_after_all_filters,1)
    y_after_all_filters_norm(k,:) = y_after_all_filters(k,:)-mean(y_after_all_filters(k,:));
    y_after_all_filters_norm(k,:) = y_after_all_filters_norm(k,:)./std(y_after_all_filters_norm(k,:));
end

[y_after_all_filters_resort_resigned,corresponding_sources,W_after_all_filters_resort] =...
    link_ICA_to_sources_and_W(y_after_all_filters_norm,sources_norm,W_after_all_filters);

figure('name','X1_after_all_filters')
nb_of_ica = size(y_after_all_filters_resort_resigned,1);
for i = 1:nb_of_ica
   subplot(nb_of_ica,2,2*i-1)
   plot(t,sources_norm(corresponding_sources(i),:))
   title(strcat('source',num2str(corresponding_sources(i))))
   subplot(nb_of_ica,2,2*i)
   plot(t,y_after_all_filters_resort_resigned(i,:))
   title(strcat('ICA signal',num2str(i)))     
end

% fourier transform of the IC related to the white noise
y_after_all_filters_resort_resigned_ft = fft(y_after_all_filters_resort_resigned(2,:));
maxft = max(abs(y_after_all_filters_resort_resigned_ft));
y_after_all_filters_resort_resignedspec = 20*log10(abs(y_after_all_filters_resort_resigned_ft)/maxft);
y_after_all_filters_resort_resigned_spec = y_after_all_filters_resort_resignedspec;

figure()
plot(freq(1:length(freq)/2),(y_after_all_filters_resort_resigned_spec(1:length(freq)/2)))
axis tight
ylabel('Magn resp [dB]')
xlabel('fequency [Hz]')
title('fequency content of the IC related to the white noise source')

%% functions

function [ICA_resort_resigned,corresponding_sources] =...
    link_ICA_to_sources(ICA_signals,sources)
% input: genormaliseerde signalen!!
% resort the ICA signals

corr_matrix = corr(ICA_signals',sources'); %rows = ica, col = sources
for m = 1:size(corr_matrix,1)
    [M,I] = max(abs(corr_matrix(m,:)));
    I_tot(m) = I;
    M_tot(m) = M;
end

new_position_in_ICA_array = 1;
for i = 1:size(sources,1)
    M_tot_copy = M_tot;
    ICA_corr_with_source = find(I_tot==i);
    M_tot_copy(ICA_corr_with_source) = 0;
    M_tot_one_source = -M_tot_copy + M_tot;
    k = length(ICA_corr_with_source);
    while k > 0
        [~,ICA_best_corr] = max(M_tot_one_source);
        ICA_resort(new_position_in_ICA_array,:) = ICA_signals(ICA_best_corr,:);
        M_tot_one_source(ICA_best_corr) = 0;
        new_position_in_ICA_array = new_position_in_ICA_array+1;
        k = k-1;
    end 
end

corr_matrix = corr(ICA_resort',sources'); %rows = ica, col = sources
% fix the sign
for m = 1:size(corr_matrix,1)
    [M,I] = max(abs(corr_matrix(m,:)));
    I_tot(m) = I;
    M_tot(m) = M;
end
for m = 1:size(corr_matrix,1)
    if corr_matrix(m,I_tot(m)) < 0
        ICA_resort(m,:) = -1*ICA_resort(m,:);
    end
end
ICA_resort_resigned = ICA_resort;
for m = 1:size(corr_matrix,1)
    [~,I] = max(abs(corr_matrix(m,:)));
    corresponding_sources(m) = I;
end
end

% same as above, but also resort the colums of the W matrix
function [ICA_resort_resigned,corresponding_sources,W_resort] =...
    link_ICA_to_sources_and_W(ICA_signals,sources,W)
% input: genormaliseerde signalen!!
% resort the ICA signals

corr_matrix = corr(ICA_signals',sources'); %rows = ica, col = sources
for m = 1:size(corr_matrix,1)
    [M,I] = max(abs(corr_matrix(m,:)));
    I_tot(m) = I;
    M_tot(m) = M;
end

new_position_in_ICA_array = 1;
for i = 1:size(sources,1)
    W_copy = W;
    M_tot_copy = M_tot;
    ICA_corr_with_source = find(I_tot==i);
    M_tot_copy(ICA_corr_with_source) = 0;
    M_tot_one_source = -M_tot_copy + M_tot;
    k = length(ICA_corr_with_source);
    while k > 0
        [~,ICA_best_corr] = max(M_tot_one_source);
        ICA_resort(new_position_in_ICA_array,:) = ICA_signals(ICA_best_corr,:);
        W_resort(new_position_in_ICA_array,:) = W_copy(ICA_best_corr,:);
        M_tot_one_source(ICA_best_corr) = 0;
        new_position_in_ICA_array = new_position_in_ICA_array+1;
        k = k-1;
    end 
end

corr_matrix = corr(ICA_resort',sources'); %rows = ica, col = sources
% fix the sign
for m = 1:size(corr_matrix,1)
    [M,I] = max(abs(corr_matrix(m,:)));
    I_tot(m) = I;
    M_tot(m) = M;
end
for m = 1:size(corr_matrix,1)
    if corr_matrix(m,I_tot(m)) < 0
        ICA_resort(m,:) = -1*ICA_resort(m,:);
    end
end
ICA_resort_resigned = ICA_resort;
for m = 1:size(corr_matrix,1)
    [~,I] = max(abs(corr_matrix(m,:)));
    corresponding_sources(m) = I;
end
end
