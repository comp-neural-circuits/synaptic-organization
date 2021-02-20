function [Sin , sCWM] = getSin(cWM , fxs , a , b)
    N = size(fxs , 3);
    sCWM = squeeze(sum(sum(fxs.*cWM , 1) , 2));
    ssCWM = a*exp(b*sCWM);
    Sin = ssCWM > rand(N,1);
end
