function [rowNCC, colNCC] = testCorr(imTemplate, imSearch)
    [M, N] = size(imTemplate);
    [R, C] = size(imSearch);

    %// Cast to double for precision
    imTemplate = im2double(imTemplate);
    imSearch = im2double(imSearch);

    neigh = im2col(imSearch, [M, N]);
    templateCol = imTemplate(:); %// Ensures we place template into single column

    imInvert = ~imTemplate;
    imInvertNoBorder = imclearborder(imInvert, 8); %// Search 8-pixel neighbourhood
    rowsToRemove = imInvertNoBorder(:) == 1;
    neigh(rowsToRemove,:) = [];

    neighMeanSubtract = bsxfun(@minus, neigh, mean(neigh));
    templateMeanSubtract = templateCol - mean(templateCol);

    numerator = bsxfun(@times, neighMeanSubtract, templateMeanSubtract);
    sumNumerator = sum(numerator);

    denominator1 = sqrt(sum(neighMeanSubtract.*neighMeanSubtract));
    denominator2 = sqrt(sum(templateMeanSubtract.*templateMeanSubtract));
    sumDenominator = denominator1 .* denominator2;

    NCC = sumNumerator ./ sumDenominator;    

    finalOutput = col2im(NCC, [M, N], [R, C]);
    finalOutput(isnan(finalOutput)) = 0;

    maxCoeff = max(abs(finalOutput(:)));
    [rowNCC, colNCC] = find(abs(finalOutput) == maxCoeff);
end