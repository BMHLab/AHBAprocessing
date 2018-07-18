% make mask for one half of matix
function half = masklHalf(A)

n = tril(nan(size(A)));
n(n==0) = 1; 
half = A.*n; 

end