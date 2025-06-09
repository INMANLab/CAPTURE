function J = calcWavTF(D,F,Fs)
% This function computes a continuous wavelet (Morlet) transform on
% a segment of EEG signal; this can be used to estimate the
%
% parameters:
% D - a row vector containing a segment of EEG signal to be transformed
% F - frequencies to sample (Hz)
% Fs - sampling rate of the time-domain signal (Hz)
% wavenumber is the size of the wavelet (typically, width=6)
%	
% returns:
% J - time-frequency transform (complex, nfreq x ntime)

wavenumber = 6;

st = 1./(2*pi*(F/wavenumber));
A = 1./sqrt(st*sqrt(pi));
J = zeros(length(F),length(D)); % initialize the time-frequency matrix
for f = 1:length(F) % loop through sampled frequencies
  t = -3.6*st(f):(1/Fs):3.6*st(f);
  m = A(f)*exp(-t.^2/(2*st(f)^2)).*exp(1i*2*pi*F(f).*t); % Morlet wavelet
  y = conv(D,m); 
  J(f,:) = y(ceil(length(m)/2):length(y)-floor(length(m)/2));
end

