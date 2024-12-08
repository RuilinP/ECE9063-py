clc;
close all;

data = edfread("RAHC7.edf");

sig1 = vertcat(data.RAHc7{:});
sig2 = vertcat(data.C131{:});

fs = 2048;

number_of_samples = length(sig1);

total_time_in_seconds = number_of_samples/fs;

number_of_10_second_intervals = floor(total_time_in_seconds/10);

number_of_samples_in_a_10_second_interval = 10*fs;

signal_over_time_T = zeros(number_of_10_second_intervals, number_of_samples_in_a_10_second_interval);
fourier_series_over_time_T = zeros(number_of_10_second_intervals, number_of_samples_in_a_10_second_interval);

%For Convenience
N = number_of_samples_in_a_10_second_interval;
amplitude_spectrum_over_time_T = zeros(number_of_10_second_intervals, N/2 + 1);

%frequency resolution in hz
frequency_resolution = fs/number_of_samples_in_a_10_second_interval;
frequency_axis = 0:frequency_resolution:fs/2;

power_spectrum_over_time_T = zeros(number_of_10_second_intervals, N/2 + 1);
delta_PSD_over_time_T = zeros(number_of_10_second_intervals,1);
theta_PSD_over_time_T = zeros(number_of_10_second_intervals,1);
alpha_PSD_over_time_T = zeros(number_of_10_second_intervals,1);
beta_PSD_over_time_T = zeros(number_of_10_second_intervals,1);


F = frequency_resolution; %again for convenience
%Note our frequency resolution F is 0.1 hz so to k hz is at k*10*F
for i = 0:number_of_10_second_intervals - 1
    index_start = i*number_of_samples_in_a_10_second_interval + 1;
    index_end = (i+1)*number_of_samples_in_a_10_second_interval;
    signal_over_time_T(i+1,:) = sig1(index_start : index_end);
    fourier_series_over_time_T(i+1,:) = fft(signal_over_time_T(i+1,:),number_of_samples_in_a_10_second_interval);
    X = fourier_series_over_time_T(i+1,:);
    A = zeros(1,N/2+1);
    A(1) = (1/N)*sqrt(real(X(1))^2 + imag(X(1))^2);
    for k = 1 : N/2 + 1
       A(k) = (2/N)*sqrt(real(X(k))^2 + imag(X(k))^2);
    end
    amplitude_spectrum_over_time_T(i+1,:) = A;
    power_spectrum_over_time_T(i+1,:) = A.^2;
    %DELTA 0.5-3.9Hz Power Spectral Density of Delta is the area under the curve between 0.5 and 4hz
    %computed using the trapezoid rule
    %note that the indexing of the frequency axis starts at 0
    delta_PSD_over_time_T(i+1,1) = trapz(frequency_axis(51:391), power_spectrum_over_time_T(i+1, 51:391));
    %THETA 4-7.9Hz ... similarly 
    theta_PSD_over_time_T(i+1,1) = trapz(frequency_axis(401:791), power_spectrum_over_time_T(i+1, 401:791));
    %Alpha 8-13.9hz 
    alpha_PSD_over_time_T(i+1,1) = trapz(frequency_axis(801:1391), power_spectrum_over_time_T(i+1, 801:1391));
    %BETA 14-37.9Hz 
    beta_PSD_over_time_T(i+1,1) = trapz(frequency_axis(1401:3791), power_spectrum_over_time_T(i+1, 1401:3791));
end
delta_theta_alpha_beta = cat(2,delta_PSD_over_time_T, theta_PSD_over_time_T, alpha_PSD_over_time_T, beta_PSD_over_time_T);
csvwrite("delta_theta_alpha_beta-RAHC7.csv", delta_theta_alpha_beta);

%Our Amplitude Spectrum A is 0 at A_k = 0
%initialize A
%subplot(4,1,1)
%A = zeros(1,N/2+1);
%A(1) = (1/N)*sqrt(real(X(1))^2 + imag(X(1))^2);
%for k = 1:N/2+1
%    A(k) = (2/N)*sqrt(real(X(k))^2 + imag(X(k))^2);
%end
plot(frequency_axis(10:300), (A(10:300).^2))
title("Amplitude Spectrum for last 10 Seconds")
xlabel("Frequency in Hz")
ylabel("Amplitude")

%first_10_seconds = signal_over_time_T(1,:);
%next_10_seconds = signal_over_time_T(2,:);
%third_10_seconds = signal_over_time_T(3,:);
%fourth_10_seconds = signal_over_time_T(4,:);



%N = length(first_10_seconds); %N is the number of samples
%Comute N=300 point fft
%X = fft(first_10_seconds, N);


