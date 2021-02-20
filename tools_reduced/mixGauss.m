function [fx] = mixGauss(X , MU , theta , S1 , S2 )

    SIGMA_P = eye(2); SIGMA_P(1,1) = S1; SIGMA_P(2,2) = S2;
    SIGMA_N = eye(2); SIGMA_N(1,1) = S1; SIGMA_N(2,2) = S2;

    ROT = [cos(theta) , -sin(theta) ; sin(theta) , cos(theta)];
    RSIGMA_P = ROT*SIGMA_P*ROT';
    RSIGMA_N = ROT*SIGMA_N*ROT';

    RMU = (ROT*[1;0])'/2;
    fx = mvnpdf(X ,MU + RMU,RSIGMA_P) - mvnpdf(X,MU  - RMU,RSIGMA_N);
    fx = fx/sum(abs(fx(:)));
end