function dauMat = getDaughterMat(tr)
    N = length(tr.D);
    dauMat = zeros(N , N);
    for ii = 1:N
        [ids , ~] = sub_tree(tr, ii);
        dauMat(ii , :) = ids;
    end   
    dauMat = dauMat + diag(ones(N,1));
end

