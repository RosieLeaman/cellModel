function equal = checkEqualityMatrices(A,B)

% check if they are the same size. Not same size, can't be equal.
if nnz(~(size(A)==size(B))) > 0
    equal = 0;
    return
end

% if they are the same size they have to have the same values

if sum(sum(A==B)) == numel(A)
    equal = 1;    
else
    equal = 0;
end