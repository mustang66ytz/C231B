%% Robust control homework
%% Verify the cnum2sys
deltaNum = complex(randn,randn) 
w = exp(3*randn)
deltaSys = cnum2sys(deltaNum,w); 
zpk(deltaSys)
[norm(deltaSys,inf) abs(deltaNum)]
[ninf,fpeak] = hinfnorm(deltaSys)
% it is observed that the infinity norm of the system is the same as the
% norm of the complex number, the infinity norm of a system is the peak
% gain of the frquency response. The system's infinity norm is at zero
% frequency, which can be proved from the code in cnum2sys.m
freqresp(deltaSys,w)