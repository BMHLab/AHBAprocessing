% make mask for one half of matix
function half = maskuHalf(A)

n = triu(nan(size(A)));
n(n==0) = 1; 
half = A.*n; 

end