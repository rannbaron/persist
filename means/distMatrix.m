function [ distMat ] = distMatrix( a, b )
%   This function creates a distance matrix for points from two persistence 
%       diagrams. The input is two matrices with occurences as rows.
%i=1;

aLen = size(a,1);
bLen = size(b,1);
%   Create matrix of pairwise distances between a and b.
distMat = pdist2(a,b);
aDist = zeros(aLen,1);
bDist = zeros(1,bLen);

%   Create vector of distances from points of matrix a to diagonal.
for i = 1:aLen
    aDist(i,1) = (a(i,2) - a(i,1))/sqrt(2);
end

%   Create block of aDist with length aLen.
aDiag = repmat(aDist,1,aLen);

%   Concatenate aDiag horizontally onto distMat.
distMat = horzcat(distMat, aDiag);

%   Create vector of distances from points of matrix b to diagonal.
for i = 1:bLen
    bDist(1,i) = (b(i,2)-b(i,1))/sqrt(2);
end

%   Add on distances for points on line y=x.
zero = zeros(1,aLen);

%   Horizontally concatenate to bDist.
bDist = horzcat(bDist,zero);

%   Create block with distance from diagonal of points on y=x and points of
%   b.
bDiag = repmat(bDist,bLen,1);

%   Finish distMat with vertical concatenation.
distMat = vertcat(distMat, bDiag);
end

