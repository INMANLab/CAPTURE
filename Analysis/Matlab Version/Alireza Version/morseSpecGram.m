function [f,coi,cfs] = morseSpecGram(Y,Fs,FPass,varargin)


tb = 30; %timebandwidth >3 and <=120 (30 mimicks morlet, 60 is default)
if nargin>3
    tb = varargin{:};
end


n = size(Y,1);
ntrials = size(Y,2);


%A gamma of 3 and timebandwidth of 30 for the morse wavelet is almost identical to the morlet
%(amor). Consider using the morse and having the ability to change timebandwidth. A higher
%timebandwidth smooths across time, which will be useful for averaging responses that are not exactly
%localized in time.
% fb = cwtfilterbank('Wavelet','amor','signallength',n,'samplingfrequency',Fs,'frequencylimits',FPass);
fb = cwtfilterbank('Wavelet','morse','signallength',n,'samplingfrequency',Fs,'frequencylimits',FPass,...
    'timebandwidth',tb);


if ~isa(Y,'double')
    Y = double(Y);
end


[cf,f,coi] = wt(fb,Y(:,1));
if ntrials<2
    cfs = cf';
    return;
end


cfs = zeros(n,length(f),ntrials); 
cfs(:,:,1) = cf';
for k=2:ntrials
    cfs(:,:,k) = wt(fb,Y(:,k))';
end
