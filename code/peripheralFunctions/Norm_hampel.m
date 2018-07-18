function dataMatrixNorm = Norm_hampel(dataMatrix)

dataMatrixNorm = zeros(size(dataMatrix));
numFeatures = size(dataMatrix,2);
        for i = 1:numFeatures % cycle through the features
            Gene = dataMatrix(:,i);
            dataMatrixNorm(:,i) = tanh_hampel(Gene);
        end

end