clear
close all
clc

addpath(genpath(pwd));
load ACC
acc=ACC(:,2);

%% Parameters (see funtion descriptions for other parameters)
Fs=2000; % sample frequency
infFreq=6; % lowest frequency analyzed (set the range for cwt and fft, discart some IMF or wavelets from wt)
supFreq=100;  % highest frequency analyzed (set the range for cwt and fft, discart some IMF or wavelets from wt)
plotFig=1; % 0 to not plot
preImpact=0.2; % the impact occured at 0.2 second after the start of the signal (400 points at 2000 Hz)
postImpact=0.25; % 0.25 seconds are considered post impact

a=100;
f=20;
c=15;
p=pi;

t=transpose(1/Fs:1/Fs:size(acc,1)/Fs);