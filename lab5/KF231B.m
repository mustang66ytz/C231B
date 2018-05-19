function [xk1k,Sxk1k,xkk,Sxkk,Sykk1,Lk,Hk,Gk,wkk] = ...
   KF231B(xkk1,Sxkk1,A,B,C,D,E,F,Swk,uk,yk)
% Implements one step of the Kalman filter, using the notation in
% slides at the end of Estimation, Part 3.  The input arguments are:
%    xkk1 = xhat_{k|k-1}
%    Sxkk1 = Sigma^x_{k|k-1}
%    A, B, C, E, F: state matrices at time k
%    Swk = variance in zero-mean disturbance w_k
%    uk = deterministic control input u_k
%    yk = measurement y_k
% The output arguments are:
%    xk1k = xhat_{k+1|k}
%    Sxk1k = Sigma^x_{k+1|k}
%    xkk = xhat_{k|k}
%    Sxkk = Sigma^x_{k|k}
%    Sykk1 = Sigma^y_{k|k-1}
%    Lk (from last slide, Part 3) for:
%           xhat_{k+1|k} = Ak x_{k|k-1} + Lk*(yk - Ck*x_{k|k-1})
%    Hk (from fact #16), for:
%           xhat_{k|k} = x_{k|k-1} + Hk*(yk - Ck*x_{k|k-1})
%    Gk (from fact #17), for: what_{k|k} = Gk*ek
%    wkk = what_{k|k}
% The input signals, xkk1, uk, yk may all be empty matrices, implying
% that the function will only update the error variances, and will not
% provide any state estimates (so xk1k, xkk and wkk will be returned
% empty as well).
% Group the variance calculations together.  Recall that the evolution of
% these quantities does not depend on measurements.

Sykk1 = C*Sxkk1*C' + F*Swk*F';
tmpK = Sxkk1*C'/Sykk1;
Sxkk = Sxkk1-tmpK*C*Sxkk1;
Gk = Swk*F'/Sykk1;
tmp = A*tmpK*F*Swk*E';
Sxk1k = A*Sxkk*A' + E*(Swk-Gk*F*Swk)*E'-tmp-tmp';

% The additional three outputs
% Notice the Lk is actually the observer gain
Hk = tmpK;
Lk = A*Hk+E*Gk;

% Calculate state estimations, if measurements and previous estimates are
% supplied. 
if ~isempty(xkk1) && ~isempty(yk)
   ek = yk - C*xkk1 - D*uk;
   xkk = xkk1 + tmpK*ek;
   xk1k = A*xkk + E*Gk*ek;
   wkk = Gk*ek;
   % If there is a control input, apply that to the evolution of the state
   % estimate.
   if ~isempty(B)
      xk1k = xk1k + B*uk;
   end
else
   % Copy empty inputs over to empty outputs
   xkk = xkk1;
   xk1k = xkk1;
   wkk = 0;
end