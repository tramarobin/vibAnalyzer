
function [intCoefficients]=interpMap(coefficients,Freq,intFreq,interpFreq,newFs,Fs)

interpTime=Fs/newFs;

[XI,YI] = meshgrid(intFreq,interpTime:interpTime:size(coefficients,2));
intCoefficients = (interp2(1:size(coefficients,2),Freq,abs(coefficients),YI,XI))';


end

