% session_5_2_1
% November 2017

clear; close all; clc; 

load ECG
dwtmode('per','nodisp'); %periodic extension for decomposition

ecg_input = ecg;

%% Wavelet decomposition
% decompose the signal using Daubechies wavelet 4 (db4)
% use functions wavedec, waverec, appcoef, detcoef
level_of_decomposition = 5;

[C,L] = wavedec(ecg_input,level_of_decomposition,'db4');
detcoef_array = detcoef(C,L,fliplr(1:level_of_decomposition));
A = appcoef(C,L,'db4',level_of_decomposition);

figure()
stem( (1:L(1)) , abs(C(1:L(1))),'.');hold on
stem( (L(1)+1:sum(L(1:2))),abs(C(L(1)+1:sum(L(1:2)))),'.' )
stem( (sum(L(1:2))+1:sum(L(1:3))),abs(C(sum(L(1:2))+1:sum(L(1:3)))),'.' )
stem( (sum(L(1:3))+1:sum(L(1:4))),abs(C(sum(L(1:3))+1:sum(L(1:4)))),'.' )
stem( (sum(L(1:4))+1:sum(L(1:5))),abs(C(sum(L(1:4))+1:sum(L(1:5)))),'.' )
stem( (sum(L(1:5))+1:sum(L(1:6))),abs(C(sum(L(1:5))+1:sum(L(1:6)))),'.' );hold off
axis tight
ylabel('magnitude coefficients')
legend('A5','D5','D4','D3','D2','D1')
title('magnitude wavelet coefficients')



%% Retain only M largest coefficient

[ecg_reconstr_10] = decomp_and_reconstr_plots(C,L,10,ecg_input,level_of_decomposition);
[ecg_reconstr_25] = decomp_and_reconstr_plots(C,L,25,ecg_input,level_of_decomposition);
[ecg_reconstr_50] = decomp_and_reconstr_plots(C,L,50,ecg_input,level_of_decomposition);
[ecg_reconstr_100] = decomp_and_reconstr_plots(C,L,100,ecg_input,level_of_decomposition);
[ecg_reconstr_500] = decomp_and_reconstr_plots(C,L,500,ecg_input,level_of_decomposition);
[ecg_reconstr_1000] = decomp_and_reconstr_plots(C,L,1000,ecg_input,level_of_decomposition);

figure()
plot((1750:2250),ecg_reconstr_50(1750:2250))
axis tight
xlabel('sample number')
title('close up of the reconstructed ecg signal M=50')


%% RMSE
[RMSE_10] = RMSE_calc(ecg_reconstr_10,ecg_input);
[RMSE_25] = RMSE_calc(ecg_reconstr_25,ecg_input);
[RMSE_50] = RMSE_calc(ecg_reconstr_50,ecg_input);
[RMSE_100] = RMSE_calc(ecg_reconstr_100,ecg_input);
[RMSE_500] = RMSE_calc(ecg_reconstr_500,ecg_input);
[RMSE_1000] = RMSE_calc(ecg_reconstr_1000,ecg_input);

M = {'M=10';'M=25';'M=50';'M=100';'M=500';'M=1000'};
RMSE = [RMSE_10;RMSE_25;RMSE_50;RMSE_100;RMSE_500;RMSE_1000];
T = table(M,RMSE)
%% Number of component


%% close up of the recontr signal and the coeff

%% Close up of the reconstructed signal

figure()
subplot(321)
plot((850:1150),ecg_reconstr_10(850:1150)); hold on
plot((850:1150),ecg_input(850:1150)); hold off
legend('reconstr','original')
axis tight
xlabel('sample number')
title('reconstr ecg, M = 10')
subplot(322)
plot((850:1150),ecg_reconstr_25(850:1150)); hold on
plot((850:1150),ecg_input(850:1150)); hold off
legend('reconstr','original')
axis tight
xlabel('sample number')
title('reconstr ecg, M = 25')
subplot(323)
plot((850:1150),ecg_reconstr_50(850:1150)); hold on
plot((850:1150),ecg_input(850:1150)); hold off
legend('reconstr','original')
axis tight
xlabel('sample number')
title('reconstr ecg, M = 50')
subplot(324)
plot((850:1150),ecg_reconstr_100(850:1150)); hold on
plot((850:1150),ecg_input(850:1150)); hold off
legend('reconstr','original')
axis tight
xlabel('sample number')
title('reconstr ecg, M = 100')
subplot(325)
plot((850:1150),ecg_reconstr_500(850:1150)); hold on
plot((850:1150),ecg_input(850:1150)); hold off
legend('reconstr','original')
axis tight
xlabel('sample number')
title('reconstr ecg, M = 500')
subplot(326)
plot((850:1150),ecg_reconstr_1000(850:1150)); hold on
plot((850:1150),ecg_input(850:1150)); hold off
legend('reconstr','original')
axis tight
xlabel('sample number')
title('reconstr ecg, M = 1000')

%% functions

function[ecg_reconstr] = decomp_and_reconstr_plots(C,L,M,original_ecg,...
    level_of_decomposition)
ecg_input = original_ecg;

% select the M largest coefficients
k = M;
C_copy = C;
largest_coeff_ind = [];
largest_coeff_value = [];
    while k > 0
        [~,I] = max(abs(C_copy));
        largest_coeff_ind = [largest_coeff_ind I];
        V = C(I);
        largest_coeff_value = [largest_coeff_value V];
        C_copy(I) = 0;
        k= k-1;
    end
    
C_new = zeros(size(C));
C_new(largest_coeff_ind) = largest_coeff_value;

% New coefficients
detcoef_array_new = detcoef(C_new,L,fliplr(1:level_of_decomposition));
A_new = appcoef(C_new,L,'db4',level_of_decomposition);

% recompose
ecg_reconstr = waverec(C_new,L,'db4');

% plots
nb_of_plots = level_of_decomposition+2;
figure('Name',['reconstruction with M = ', num2str(M)])
subplot(nb_of_plots,1,1)
plot(ecg_input); hold on
plot(ecg_reconstr);hold off
axis tight
title(['comparison original ecg and reconstruction, M=', num2str(M)])
xlabel('sample number')
legend('original','reconstructed')
subplot(nb_of_plots,1,2)
indx_A = find(A_new); % only plot the non-zero elements
stem(indx_A,A_new(indx_A),'.');% only plot the non-zero elements
axis tight
ylabel(['A', num2str(level_of_decomposition)])
% title(['appcoef #', num2str(level_of_decomposition)])
set(gca,'XTickLabel',[])
for k = 1:level_of_decomposition
    subplot(nb_of_plots,1,k+2);
    indx = find(cell2mat(detcoef_array_new(k)));% only plot the non-zero elements
    array = cell2mat(detcoef_array_new(k));% only plot the non-zero elements
    stem(indx,array(indx),'.');% only plot the non-zero elements
    axis tight
    ylabel(['D', num2str(level_of_decomposition-k+1)])
%     title(['detcoef #', num2str(level_of_decomposition-k+1)])
    set(gca,'XTickLabel',[])
end
end


function[] = decomp_plot(C,L,ecg_input,level_of_decomposition)

detcoef_array = detcoef(C,L,fliplr(1:level_of_decomposition));
A = appcoef(C,L,'db4',level_of_decomposition);

% plots
nb_of_plots = level_of_decomposition+2;
figure('Name',['wavelet decomposition of level', num2str(level_of_decomposition)])
subplot(nb_of_plots,1,1)
plot(ecg_input);
axis tight
title('original ecg')
xlabel('sample number')
subplot(nb_of_plots,1,2)
indx_A = find(A); % only plot the non-zero elements
stem(indx_A,A(indx_A),'.');% only plot the non-zero elements
axis tight
ylabel(['A', num2str(level_of_decomposition)])
% title(['appcoef #', num2str(level_of_decomposition)])
set(gca,'XTickLabel',[])
for k = 1:level_of_decomposition
    subplot(nb_of_plots,1,k+2);
    indx = find(cell2mat(detcoef_array(k)));% only plot the non-zero elements
    array = cell2mat(detcoef_array(k));% only plot the non-zero elements
    stem(indx,array(indx),'.');% only plot the non-zero elements
    axis tight
    ylabel(['D', num2str(level_of_decomposition-k+1)])
%     title(['detcoef #', num2str(level_of_decomposition-k+1)])
    set(gca,'XTickLabel',[])
end
end


function[] = decomp_plot_2(C,L,ecg_input,level_of_decomposition)

detcoef_array = detcoef(C,L,fliplr(1:level_of_decomposition));
A = appcoef(C,L,'db4',level_of_decomposition);


% plots
nb_of_plots = level_of_decomposition+2;
figure('Name',['wavelet decomposition of level', num2str(level_of_decomposition)])
subplot(nb_of_plots,1,1)
plot(ecg_input);
axis tight
title('original ecg')
xlabel('sample number')
subplot(nb_of_plots,1,2)
indx_A = find(A); % only plot the non-zero elements
plot(indx_A,A(indx_A));% only plot the non-zero elements
axis tight
ylabel(['A', num2str(level_of_decomposition)])
% title(['appcoef #', num2str(level_of_decomposition)])
set(gca,'XTickLabel',[])
for k = 1:level_of_decomposition
    subplot(nb_of_plots,1,k+2);
    indx = find(cell2mat(detcoef_array(k)));% only plot the non-zero elements
    array = cell2mat(detcoef_array(k));% only plot the non-zero elements
    plot(indx,array(indx));% only plot the non-zero elements
    axis tight
    ylabel(['D', num2str(level_of_decomposition-k+1)])
%     title(['detcoef #', num2str(level_of_decomposition-k+1)])
    set(gca,'XTickLabel',[])
end
end

function [RMSE] = RMSE_calc(signal,reference)
RMSE = sqrt(mean((signal - reference).^2));  % Root Mean Squared Error
end
