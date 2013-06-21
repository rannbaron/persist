function [ mean ] = kPersistenceMean( cell_diagrams, epsilon, zigzag, lp )
% Given a cell-array in which each cell is an n-by-2 persistence diagram,
% this program returns a k-by-2 array with the mean of the given diagrams.
%	cell_diagrams: a 1-by-n cell array, each of which contains a p-by
%		2 persistence diagram.
%	epsilon: the threshold. Points with persistence (y-x) less than
%		epsilon will be ignored.
%   zigzag: the "margin for error" when checking for equality /
%       convergence. If Y and Y1 have the same size and the arrays differ
%       pointwise by no more than 'zigzag', then they will be considered
%       equal and the algorithm will terminate.

for i=1:size(cell_diagrams, 2)                                  % First, we iterate through the provided diagrams, replacing empty diagrams --
    if size(cell_diagrams{i},1)==0&&size(cell_diagrams{i},2)==0 % which are by-default represented as 0-by-0 arrays --
        cell_diagrams{i}=zeros(0,2);                            % instead with 0-by-2 arrays, which are bounds-compatible later on.
    end
end
for i=1:size(cell_diagrams,2)                                       % This loop checks for, and responds accordingly to, points either below the diagonal or below the threshold.
    k=size(cell_diagrams{i},1);
    j=1;
    while j<=k                                                      % A while loop is used instead of a for loop: Matlab dislikes concurrent modification of array size within loops.
        if cell_diagrams{i}(j,2)-cell_diagrams{i}(j,1)<0            % Below the diagonal.
            err = MException('ResultChk:OutOfRange', 'Invalid diagram(s). Cannot have points below the diagonal');
            throw(err)
        elseif cell_diagrams{i}(j,2)-cell_diagrams{i}(j,1)<epsilon 	% Below the persistency threshold.
            cell_diagrams{i}(j,:)=[];                               % The row is deleted.
            j=j-1;                                                  % Current index must decrease, as well as
            k=k-1;                                                  % the overall loop size.
        end
        j=j+1;
    end
end

% The algorithm is a "fixed-point" algorithm, iterating until the
% diagram 'Y' stabilizes. The outermost loop uses a "conveyor belt"
% strategy: a "new" array 'Y1' is compared to the "old" array Y; if changes
% have taken place, the new array retreats to the position of the old,
% undergoes further changes, and is compared once again. When the two
% arrays are identical, the algorithm terminates.

G=zeros(0);             % 'G' will store "best bijections" between 'Y' and each diagram 'X_j'. Rows of 'G' will corresond to matchings.

Y=[-1, -2] ;            % The old array is initialized to a dummy: purposefully different.
Y1=cell_diagrams{1};    % The new array is initialized to the first provided persistence diagram.
while ~equals(Y, Y1, zigzag)    % Main check of equality for fixed-point algorithm. Calls the "equals" subroutine for equality which supports floating-point discrepancies, and passes in the "standard for equality" 'num'.
    Y=Y1;                       % First step in the loop: the "new" becomes the "old" as we prepare for a further round of changes.
    Y1=zeros(0,2);              % 'Y1' is now emptied; it will be filled as the loop progresses.
    
    k=0; % 'k' will later be used in indexing to compensate for empty rows.
    for j=1:size(cell_diagrams, 2)                  % Now, we fill 'G' with best bijections. The jth column of 'G' will denote the best bijection between 'Y' and 'X_j'.
        %dist = distMatrix(Y, cell_diagrams{j}).^k;     % Pairwise distances between 'Y' and 'X_j'. The subroutine 'DistMatrix' supports differences in array sizes by padding with matches to the diagonal.
        %[m c] = m_hungarian(dist);                  % The subroutine 'm_hungarian' uses the distance matrix to produce a best bijection between its rows (Y) and its columns (X_j).
        [m c] = kWass(Y, cell_diagrams{j},lp);
        %c=c^(1/k);
        m=m+1;                                      % Computer scientists count from 0; mathematicians count from 1.
        G(1:size(m,1),j)=m;                         % The hungarian matching from 'Y' to 'X_j' is copied to the appropriate column of 'G'.
    end
    
    extras = zeros(0,2);    % 'extras' will be used to store matchings from the diagonal to non-diagonals. These points are "drawn out of the diagonal" and "come from nowhere"; they'll later be concatenated to Y1.
    e=1;                    % 'e' will be used to index the array 'extras'.
    
    for i=1:size(G,1)   % This large, important loop crawls down the rows of 'G', computing, for each row, an arithmetic mean of its elements and storing that mean into the appropriate index of either 'Y' or 'extras'.
        point_count=0;                                          % 'point_count' is used in the small loop below to identify a row consisting entirely of diagonals.
        for j=1:size(G,2)                                       % Does the ith row contain non-diagonal points? We determine whether a point is "actual" by comparing the destination number to the size of the original diagram.
            if G(i,j)<=size(cell_diagrams{j},1)&&G(i,j)~=0      % "In-bounds" status is determined by comparing the value of G(i,j) to the total size of the original cell diagram X_j to which the matching refers.
                point_count=1;                                  % Denote the discovery of a "healthy" in-bounds point
                break;                                          % and break from the loop.
            end
        end
        if point_count==0                                       % If no in-bounds points have been discovered, we've encountered a row of matches entirely to the diagonal.
            k=k+1;                                              % We must mark the "lost point" through the indexing variable, 'k'. The necessity of this will become evident later.
            continue;                                           % Having documented the lost point, we may skip the rest of the loop and continue to the next row of 'G'.
        end
        temp=zeros(1,2);    % The variable 'temp' will store the sum, and eventually the arithmentic mean, of the matching indicated by the ith row of 'G'.
        if i<=size(Y,1)         % FROM: Are we matching from a real point, or from a diagonal? If 'i' is within the range of 'Y', we're matching from an ACTUAL point.    
            for j=1:size(G,2)
                if G(i,j)<=size(cell_diagrams{j},1)     % TO: Are we matching to a real point, or to the diagonal? If G(i,j) is within the range of the diagram 'X_j', we're matching to an ACTUAL point.
                    temp=temp+cell_diagrams{j}(G(i,j),:);       % Using the value of G(i,j), pull the appropriate off-diagonal point from the original persistence diagram 'X_j' and add it to 'temp'.
                else                                    % TO: If G(i,j) is outside of the range of 'X_j', then we're matching to a DIAGONAL point.
                    off=Y(i,:);                                             % To calculate the position of the desired "destination" diagonal point, first, we retrieve the position of the "source" point from 'Y1',
                    diag=[(off(1,1)+off(1,2))/2 (off(1,1)+off(1,2))/2];     % and then, we "average the coordinates" to project it onto the diagonal.
                    temp=temp+diag;                                         % Finally, we add the diagonal destination to the running total 'temp'.
                end
            end
            temp=temp*(1/size(G,2));                    % We convert 'temp' from a sum to an arithmetic mean by dividing by the number of columns of 'G'.
            Y1(i-k,:)=temp;                             % Finally, we store 'temp' into 'Y1'. Because empty rows correspond to "lost points", we must subtract the amount of lost points, 'k', from our index 'i'.
        else                    % FROM: If 'i' is outside of the range of 'Y', we're matching from a DIAGONAL point.
            point_count=0;          % 'point_count' will be used, in the loop below, to determine how many of the points in the ith row are "actual" and not diagonal.
            for j=1:size(G,2)   % Thus, below, we have understood matching from diagonals. The non-diagonals in the matching are averaged, the average is projected to the diagonal, and finally the projection is used to represent the diagonals.
                if G(i,j)<=size(cell_diagrams{j},1)&&G(i,j)~=0              % "In-bounds-ness" is calculated, as above.
                    temp=temp+cell_diagrams{j}(G(i,j),:);                   % If an in-bounds point is discovered, its location is added to the total temp,
                    point_count=point_count+1;                              % and 'point_count' is incremented.
                end                                                         % Note: we know that 'point_count' could not be zero here; if it were, we would already have "continue"-ed on line 51.
            end
            off=temp/point_count;                                           % The sum of off-diagonal points is averaged,
            diag=[(off(1,1)+off(1,2))/2 (off(1,1)+off(1,2))/2];             % and then projected to the diagonal.
            
            temp=zeros(1,2);                                                % 'temp' is emptied.
            for j=1:size(G,2)                                               % Finally, we loop through again, using are "representative" to stand for diagonal points and treating the rest naturally.
                if G(i,j)<=size(cell_diagrams{j},1)                         % Check for in-bounds-ness.
                    if G(i,j)~=0                                            % A G(i,j) value of zero signifies a meaningless "padded" value outside of the range even of the hungarian matching between 'Y' and 'X_j'.
                        temp=temp+cell_diagrams{j}(G(i,j),:);               % Add a non-diagonal point.
                    end
                else
                    temp=temp+diag;                                         % Add a diagonal point, using its representative.
                end
            end
            temp=temp*(1/size(G,2));                                        % Finally, the sum is turned into an arithmetic mean.
            extras(e,:)=temp;                                               % The mean is stored into the 'extras' array,
            e=e+1;                                                          % and the indexing variable 'e' is incremented.
        end
    end
    Y1=vertcat(Y1,extras);      % The grand finale: the now-filled 'Y1' array is concatenated with the new arrivals in the 'extras' array. 'Y1' is now ready for another comparison.
    
    q = size(Y1,1); % Now, we scan through the completed Y1 and remove, according to the threshold value, points with inadequate persistence.
    i=1;
    while i<=q
        if Y1(i,2)-Y1(i,1)<epsilon
            Y1(i,:)=[];
            i=i-1;
            q=q-1;
        end
        i=i+1;
    end
end
mean=Y1; % After termination of the giant while loop, the return result is stored.
end

function [ answer ] = equals(Y, Y1, zigzag)
%This subroutine checks for equality.
answer = true;
if size(Y,1)~=size(Y1,1)
    answer=false;
    return;
end

for i=1:size(Y,1)
    if norm(Y(i,:)-Y1(i,:))>zigzag
        if abs(Y(i,1)-Y(i,2))>zigzag||abs(Y1(i,1)-Y1(i,2))>zigzag
            answer = false;
        end
    end
end
end