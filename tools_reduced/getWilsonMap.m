function cmap = getWilsonMap(N)
    rawMap = [213, 71, 146; 106, 84, 158; 64, 93, 166; 115, 199, 217; 116, 183, 72; 147, 193, 68; 248, 189, 22; 229, 42, 43; 213, 71, 146]/256;
    M = size(rawMap , 1);
    cmap = zeros(N , 3);
    cmap(: , 1) = interp1((0:M-1)/(M-1) , rawMap(:,1) , (0:N-1)/(N-1));
    cmap(: , 2) = interp1((0:M-1)/(M-1) , rawMap(:,2) , (0:N-1)/(N-1));
    cmap(: , 3) = interp1((0:M-1)/(M-1) , rawMap(:,3) , (0:N-1)/(N-1));

end