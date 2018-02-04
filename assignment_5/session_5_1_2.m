% session_5_1_2.m
% November 2017

clear; close all; clc; 
addpath ./FastICA_25/

%% Data
load EEG.mat; 
fs = 128;
t= (0:1/fs:20-1/fs);
freq = linspace(0,fs,length(t));


% plot the EEG data using eegplot_simple
figure()
eegplot_simple(data,channels)
title('EEG original multichannel signal')
%% Apply fastICA

[ICA_data,A,W] = fastica(data); % dit maar 1 keer runnen want kan elke run anders zijn





%% Artefact removal 
% based on the temporal and spatial information obtained by ICA, select which components to remove. 
% Adjust the mixing matrix accordingly and use it to construct a version of the EEG without artifacts. 
% Plot the clean EEG  compare to the first plot in this exercise.

% plot the temporal information using eegplot_s figure()
for i = 1:size(ICA_data,1)
    subplot(size(ICA_data,1),1,i);
    plot(t, ICA_data(i,:));
    title(['icasig signal #', num2str(i)]);
    xlabel('Time (s)');
end


% % plot the topographic information using topoPlot
figure()
topoPlot(A) % moet hier A of W


% find the blink channel (dit is valsspelen, maar das gwn zodat de code
% altijd runt)
load('blink_whole_channel') % opgeslagen uit een vorige run
corr_matrix = corr(ICA_data',blink_channel');

[~,IC_blink_position] = max(abs(corr_matrix)); %deze component moet verwijderd worden



%% remove the selected component(s)

% de te verwijderen kolom kan verschillen van run tot run!


W_inv = pinv(W); % pseudoinverse cause W does not has to be square
W_inv(:,IC_blink_position) = 0;
X_without_blink = W_inv*ICA_data;


figure()
eegplot_simple(X_without_blink,channels)
title('EEG original multichannel eyeblink removed')



[ICA_data_without_blink,~,W_whitout_blink] = fastica(X_without_blink);
for k = 1:size(ICA_data_without_blink,1)
    ICA_data_without_blink(k,:) = ICA_data_without_blink(k,:)-mean(ICA_data_without_blink(k,:));
    ICA_data_without_blink(k,:) = ICA_data_without_blink(k,:)./std(ICA_data_without_blink(k,:));
end

figure()
for i = 1:size(ICA_data_without_blink,1)
    subplot(size(ICA_data_without_blink,1),1,i);
    plot(t,ICA_data_without_blink(i,:));
    title(['icasig signal blink removed #', num2str(i)]);
    xlabel('Time (s)');
end


%% Functions

function [CCF] = calc_CCF(template_input,signal_input)
M = length(template_input);
N = length(signal_input);
signal_input_add = [signal_input zeros(1,M)]; %nullen toevoegen zodat de CCF even lang is als signal
CCF = [];
for k = 1:N % omdat er nullen achteraan signal staan, (als ni: N-M)
    cc = 0;
    for m = 1:M
        cc = cc + template_input(m)*signal_input_add(m+k);
    end    
    CCF = [CCF cc];
end
end