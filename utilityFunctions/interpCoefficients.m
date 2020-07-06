function intCoefficients=interpCoefficients(coefficients,f,intFreq,interpFreq,newFs,Fs,indIntInfLim,indIntSupLim,newPreImpactPoints,newPostImpactPoints)


[intCoefficients]=interpMap(coefficients,f,intFreq,interpFreq,newFs,Fs);
    
indices=min([indIntInfLim indIntSupLim]):max([indIntInfLim indIntSupLim]);

intCoefficients=intCoefficients(indices,round(newPreImpactPoints:newPreImpactPoints+newPostImpactPoints-1));


end

