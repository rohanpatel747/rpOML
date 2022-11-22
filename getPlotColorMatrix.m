function colorIndex = getPlotColorMatrix(colorType, A, sortbyIndex)
%GETPLOTCOLORMATRIX Returns RGB Triplet of Intensity in Row of A
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. A          [nx9 Dbl.]  Col. Matrix of Data
%       2. sortbyIndex[1x1 Dbl.]  Row to Color Data to Sort By:
%       3. colorType  [str.]      MATLAB Colormap Type:
%                                   'parula','cool','jet', etc.
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output:
%       1. colordata  [nx1  Dbl.] RGB color pair w.r.t A(:,sortbyIndex)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
    minval = min(A(:,sortbyIndex));
    %maxval = max(A(:,sortbyIndex));
    
    B = A(:,sortbyIndex);
    B = B-minval;
    B = (1/max(B)).*B;
    
    cmap     = colormap(colorType);
    cmap     = [((1:length(cmap)).'./length(cmap)) , cmap];
    ncolors  = length(cmap);
    clrIdx   = ((1/ncolors).*(1:ncolors)).';
    
    d        = abs(bsxfun(@minus, B(:).', clrIdx(:))); 
    [~, ind] = min(d);
    clrIdxA  = clrIdx(ind);
    
    A = [clrIdxA, A];
    for i=1:length(cmap)
    
        ci = cmap(i);
        ai = find(A(:,1)==ci);
        for j=1:length(ai)
            ind = ai(j);
            colorIndex(ind,1:3) = cmap(i,2:4);
        end
    end
end