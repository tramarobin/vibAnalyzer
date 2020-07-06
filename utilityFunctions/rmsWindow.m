%% OUPUT
% Signal rms

%% INPUTS
% Signal (time=column)
% Window (seconds)
% Hz

function Signal_rms=rmsWindow(Signal,Window,Fs)

Points=Window*Fs;
Signal_rms=sqrt(movmean(Signal.^2,[Points/2 Points/2]));


end