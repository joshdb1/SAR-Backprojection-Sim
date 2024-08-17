% Computes a chirp signal with parameters requested
% Written by Josh Brock 2021
%
% Inputs:
%   t  : the time (or times) at which the function is evaluated
%   T  : the pulse duration
%   K  : the linear FM rate
%   fc : the center frequency
%
% Outputs:
%   chirp : the output chirp signal. Length of chirp matches length of t

function chirp = mychirpcausal(t,T,K,fc)
    chirp = rectangularPulse(0,T,t).*exp(1i*(2*pi*fc.*(t-T/2)+pi*K*(t-T/2).^2));
end