% make mask for one half of matix
function half = maskHalf(A)

n = triu(nan(size(A)));
n(n==0) = 1; 
half = A.*n; 

end