function modes=peemd(x,Nstd,NR,MaxIter, Hlim, n, tau)

% WARNING: for this code works it is necessary to include in the same
%directoy the file emd.m developed by Rilling and Flandrin.
%This file is available at %http://perso.ens-lyon.fr/patrick.flandrin/emd.html
%We use the default stopping criterion.
%We use the last modification: 3.2007
%
% This version was run on Matlab 7.10.0 (R2010a)
%----------------------------------------------------------------------
%   INPUTs
%   x: signal to decompose
%   Nstd: noise standard deviation
%   NR: number of realizations
%   MaxIter: maximum number of sifting iterations allowed.
%   Hlim : minimum permutation enthropy to continue CEEMD
%
%  OUTPUTs
%  modes: contain the obtained modes in a matrix with the rows being the modes
%   its: contain the sifting iterations needed for each mode for each realization (one row for each realization)
% -------------------------------------------------------------------------
%  Syntax
%
%  modes=ceemdan(x,Nstd,NR,MaxIter)
%  [modes its]=ceemdan(x,Nstd,NR,MaxIter)
%
%--------------------------------------------------------------------------
% This algorithm was presented at ICASSP 2011, Prague, Czech Republic
% Plese, if you use this code in your work, please cite the paper where the
% algorithm was first presented.
% If you use this code, please cite:
%
% M.E.TORRES, M.A. COLOMINAS, G. SCHLOTTHAUER, P. FLANDRIN,
%  "A complete Ensemble Empirical Mode decomposition with adaptive noise,"
%  IEEE Int. Conf. on Acoust., Speech and Signal Proc. ICASSP-11, pp. 4144-4147, Prague (CZ)
%
% -------------------------------------------------------------------------
% Date: June 06,2011
% Authors:  Torres ME, Colominas MA, Schlotthauer G, Flandrin P.
% For problems with the code, please contact the authors:
% To:  macolominas(AT)bioingenieria.edu.ar
% CC:  metorres(AT)santafe-conicet.gov.ar
% -------------------------------------------------------------------------

% Default values
if nargin < 7, tau = 1; end
if nargin < 6, n = 3; end
if nargin < 5, Hlim = 1.8; end

x=x(:)';
desvio_x=std(x);
x=x/desvio_x;
NR = floor(NR/2);
modes=zeros(size(x));
temp=zeros(size(x));
aux=zeros(size(x));
acum=zeros(size(x));
iter=zeros(NR,round(log2(length(x))+5));

for i=1:NR
    white_noise{i}=randn(size(x));%creates the noise realizations
end

for i=1:NR
    modes_white_noise{i}=emdCustom(white_noise{i});%calculates the modes of white gaussian noise
end

for i=1:NR %calculates the first mode
    temp_p=x+Nstd*white_noise{i};
    temp_p=emdCustom(temp_p,'MAXMODES',1,'MAXITERATIONS',MaxIter);
    temp_p=temp_p(1,:);
    aux=aux+temp_p/(2*NR);
    temp_m=x-Nstd*white_noise{i};
    temp_m=emdCustom(temp_m,'MAXMODES',1,'MAXITERATIONS',MaxIter);
    temp_m=temp_m(1,:);
    aux=aux+temp_m/(2*NR);
end

modes=aux; %saves the first mode
k=1;
aux=zeros(size(x));
acum=sum(modes,1);
H = pentropy(modes,n,tau); %permutaion entropy of the first mode
while  nnz(diff(sign(diff(x-acum))))>2 && H>Hlim %calculates the rest of the modes
    for i=1:NR
        tamanio=size(modes_white_noise{i});
        if tamanio(1)>=k+1
            noise=modes_white_noise{i}(k,:);
            noise=noise/std(noise);
            noise=Nstd*noise;
            try
                [temp_p, o, it]=emdCustom(x-acum+std(x-acum)*noise,'MAXMODES',1,'MAXITERATIONS',MaxIter);
                temp_p=temp_p(1,:);
            catch
                it=0;
                temp_p=x-acum;
            end
            try
                [temp_m, o, it]=emdCustom(x-acum-std(x-acum)*noise,'MaxNumIMF',1,'MAXITER',MaxIter);
                temp_m=temp_m(1,:);
            catch
                it=0;
                temp_m=x-acum;
            end
            temp = (temp_p+temp_m)/2;
        else
            [temp, o, it]=emdCustom(x-acum,'MaxNumIMF',1,'MAXITER',MaxIter);
            temp=temp(1,:);
        end
        aux=aux+temp/NR;
    end
    H = pentropy(aux,n,tau);
    modes=[modes;aux];
    aux=zeros(size(x));
    acum=zeros(size(x));
    acum=sum(modes,1);
    k=k+1;
end
[temp]=emdCustom(x-acum,'MAXITERATIONS',MaxIter);
nb_mode_emd=length(temp(:,1));
% iter=[iter it];
modes=[modes;temp];
[a b]=size(modes);
iter=iter(:,1:a);
modes=modes*desvio_x;
its=iter;