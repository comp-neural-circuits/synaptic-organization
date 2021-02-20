function [RGL , LB , UB , DX , CCSAMP] = bootstrapRegression(X , Y , N)
    if nargin < 3
        N = 1000;
    end
    DX = linspace(min(X) , max(X) , 100);
    f = fit(X , Y , 'poly1');
    RGL = f(DX);
    REGLSAMP = zeros(N , length(DX));
    CCSAMP = zeros(N,1);
    for xx = 1:N
        BTSTR_ID = datasample(1:length(X) , length(X));
        BTSTR_X = X(BTSTR_ID); BTSTR_Y = Y(BTSTR_ID);
        BTSTR_f = fit(BTSTR_X , BTSTR_Y , 'poly1');
        REGLSAMP(xx , :) = BTSTR_f(DX);
        tmp = corrcoef(BTSTR_X , BTSTR_Y);
        CCSAMP(xx) = tmp(1,2);
    end
    LB = prctile(REGLSAMP , 5); UB = prctile(REGLSAMP , 95);
end