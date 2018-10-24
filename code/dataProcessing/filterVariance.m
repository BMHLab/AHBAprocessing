function INDkeep = filterVariance(DataTable, scale, percFilter)
% variance will be different across subjects and high variance on
% non-normalised combine data can result purely from differences between
% subjects, therefore we calculate variance for each subject and remove
% probes that are consistently low variance
varVALS = cell(size(DataTable,1),1);
nrGenes = round(size(DataTable.Expression{1},1)*percFilter/100);
LOWvar = zeros(size(DataTable.Expression{1},1), size(DataTable,1));

for s=1:size(DataTable,1)
    switch scale
        case 'log2'
            expVals = DataTable.Expression{s};
        case 'normal'
            expVals = 2.^(DataTable.Expression{s});
    end
    varSUBJ = var(expVals,0,2);
    varVALS{s} = varSUBJ;
        
    [ ~, ix ] = sort(varSUBJ, 'ascend' );
    
    for ii=1:size(expVals,1)
        LOWvar(ii,s) = ix(ii);
    end
end
% vectorise those genes that are kept
A = LOWvar(1:nrGenes,:); 
A = A(:); 
% find indexes that are repeated in all 6 cubjects and save them for
% removal
[n] = histc(A, unique(A));
multiple = n == 6;
index    = unique(A(multiple));
% get the indexes of probes to keep (the opposite of index)
INDkeep = setdiff(1:size(DataTable.Expression{1},1), index)'; 

end