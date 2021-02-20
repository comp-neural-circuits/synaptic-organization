function rdMat = augDMat(dMat , dauMat , subpos)
    N = length(subpos);
    rdMat = zeros(N , N);
    for ii = 1:N
        for jj = (ii+1):N
             if dauMat(ii , jj) == 2
                rdMat(ii , jj) = abs(subpos(ii) - subpos(jj));
             elseif dauMat(ii , jj) == 1
                rdMat(ii , jj) = dMat(ii , jj) - subpos(ii) + subpos(jj);
             else
                rdMat(ii , jj) = dMat(ii , jj) + subpos(ii) - subpos(jj);
             end
        end
    end
    rdMat = rdMat + rdMat'+ diag(ones(N,1));
end

