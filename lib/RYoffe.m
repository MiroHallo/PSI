function f = RYoffe(N,dt,t0,tR,tP,D)
%---------------------------------------------------------------------
%  Generates discrete time array with Regularized Yoffe function
%  (Result in RYoffe is slip-rate time function)
%  Following: Tinti E., Fukuyama E., Piatanesi A., Cocco M. (2005). A Kinematic Source-Time
%             Function Compatible with Earthquake Dynamics, BSSA, 95(4), 1211–1223.
%
%  Code author: Miroslav Hallo
%  Charles University in Prague, Faculty of Mathematics and Physics
%  E-mail: hallo@karel.troja.mff.cuni.cz
%  Revision 1/2018: The first version of the function
%
%  Copyright (C) 2018  Miroslav Hallo
%
%  This program is published under the GNU General Public License (GNU GPL).
%  
%  INPUT: N .. number of samples
%         dt .. sampling interval [sec]
%         t0 .. rupture time [sec] when the 1st sample is 0.0 time
%         tR .. rise time [sec] relatively to t0
%         tP .. peak time [sec] relatively to t0
%         D .. total slip [m]
%---------------------------------------------------------------------

f = zeros(N,1);

%----------------------------------
%  Prepare tS
tS = tP/1.27;
if(tR<2*tS)
    f(1)=-12345;
    return
end

%----------------------------------
%  Equation A13
K = D*2/(pi*tR*tS^2);
tfc = dt*[0:N-1]-t0;

for i=1:N
    t = tfc(i);
    if(t<0)
        continue
    elseif(t<tS)
        C1=(0.5*t+0.25*tR)*sqrt(t*(tR-t))+(t*tR-tR^2)*asin(sqrt(t/tR))-(3/4*tR^2)*atan(sqrt((tR-t)/t));
        C2=(3/8)*pi*tR^2;
        f(i)=K*(C1+C2);
    elseif(t<2*tS)
        C1=(0.5*t+0.25*tR)*sqrt(t*(tR-t))+(t*tR-tR^2)*asin(sqrt(t/tR))-(3/4*tR^2)*atan(sqrt((tR-t)/t));
        C2=(3/8)*pi*tR^2;
        C3=(tS-t-0.5*tR)*sqrt((t-tS)*(tR-t+tS))+tR*2*(tR-t+tS)*asin(sqrt((t-tS)/tR))+(3/2)*(tR^2)*atan(sqrt((tR-t+tS)/(t-tS)));
        f(i)=K*(C1-C2+C3);
    elseif(t<tR)
        C1=(0.5*t+0.25*tR)*sqrt(t*(tR-t))+(t*tR-tR^2)*asin(sqrt(t/tR))-(3/4*tR^2)*atan(sqrt((tR-t)/t));
        C3=(tS-t-0.5*tR)*sqrt((t-tS)*(tR-t+tS))+tR*2*(tR-t+tS)*asin(sqrt((t-tS)/tR))+(3/2)*(tR^2)*atan(sqrt((tR-t+tS)/(t-tS)));
        C4=(-tS+0.5*t+0.25*tR)*sqrt((t-2*tS)*(tR-t+2*tS))-tR*(tR-t+2*tS)*asin(sqrt((t-2*tS)/tR))-(3/4)*(tR^2)*atan(sqrt((tR-t+2*tS)/(t-2*tS)));
        f(i)=K*(C1+C3+C4);
    elseif(t<tR+tS)
        C3=(tS-t-0.5*tR)*sqrt((t-tS)*(tR-t+tS))+tR*2*(tR-t+tS)*asin(sqrt((t-tS)/tR))+(3/2)*(tR^2)*atan(sqrt((tR-t+tS)/(t-tS)));
        C4=(-tS+0.5*t+0.25*tR)*sqrt((t-2*tS)*(tR-t+2*tS))-tR*(tR-t+2*tS)*asin(sqrt((t-2*tS)/tR))-(3/4)*(tR^2)*atan(sqrt((tR-t+2*tS)/(t-2*tS)));
        C5=(pi/2)*tR*(t-tR);
        f(i)=K*(C5+C3+C4);
    elseif(t<tR+2*tS)
        C4=(-tS+0.5*t+0.25*tR)*sqrt((t-2*tS)*(tR-t+2*tS))-tR*(tR-t+2*tS)*asin(sqrt((t-2*tS)/tR))-(3/4)*(tR^2)*atan(sqrt((tR-t+2*tS)/(t-2*tS)));
        C6=(pi/2)*tR*(2*tS-t+tR);
        f(i)=K*(C4+C6);
    else
        break
    end
end

return



