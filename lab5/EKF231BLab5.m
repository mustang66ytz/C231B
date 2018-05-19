function [xk1k,Sxk1k,xkk,Sxkk,Sykk1,Lk,Hk,Gk,wkk] = ...
   EKF231BLab5(xkk1,Sxkk1,f,h,A,B,C,E,F,Swk,uk,yk)

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
   ek = yk -h(xkk1, uk);
   xkk = xkk1 + tmpK*ek;
   xk1k = f(xkk, uk) + E*Gk*ek;
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
