%% INPUTS
% Coefficients
% Frequencies
% Sample Frequency
% New figure ?

%% OUTPUT
% Figure of the wavelet transform

function []=displayCoeffMap(Coeff,F,Time,varargin)
%% Optional inputs
p = inputParser;

% utilities
addParameter(p,'newFig',1,@isnumeric); % 1 = new figure
addParameter(p,'nXticks',6,@isnumeric); % ticks for x
addParameter(p,'nYticks',6,@isnumeric); % ticks for y
addParameter(p,'colorbarOk',0,@isnumeric); % 1 to plot colorbar
addParameter(p,'units','m\cdots^-^2',@ischar); % units
addParameter(p,'titleFig',[],@ischar); % units

parse(p,varargin{:});

newFig=p.Results.newFig;
nx=p.Results.nXticks;
ny=p.Results.nYticks;
colorbarOk=p.Results.colorbarOk;
units=p.Results.units;
titleFig=p.Results.titleFig;

if newFig==1
    figure
end

imagesc(Time,F,flipud(abs(Coeff)))
ylabel('Frequency (Hz)')
colormap(cbrewer('seq','Reds', 64))
yticks(linspace(min(F),max(F),ny))
if max(F)<1
    ylab=linspace(min(F),max(F),ny);
    for i=1:ny
        ylabs{i}=sprintf('%0.3f',ylab(i));
    end
    yticklabels(fliplr(ylabs))
    
elseif max(F)>20
    yticklabels(round(linspace(max(F),min(F),ny)))
else
    yticklabels(linspace(max(F),min(F),ny))
end
xticks(linspace(min(Time),max(Time),nx))
xlab=linspace(0,max(Time),nx);
if max(Time)<10
    for i=1:nx
        xlabs{i}=sprintf('%0.2g',xlab(i));
    end
else
    for i=1:nx
        xlabs{i}=sprintf('%0.0f',xlab(i));
    end
end
xticklabels(xlabs)
set(gca,'Box','off');
xlabel('Time (s)');

if colorbarOk
    Co=colorbar('EastOutside');
    Co.Label.String=['Amplitude (' units ')'];
end
title(titleFig)
end
