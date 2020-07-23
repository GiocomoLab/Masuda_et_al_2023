function ind = sortByTwoCols(AStruct)

% Input:  the 'allCells' struct
% Return: sorting index 'ind',
%         which is a size(allCells.peakiness,1)-by-1 vector

Col1 = nanmean(AStruct.peakiness,2);
Col2 = nanmean(mean(AStruct.spatialFR10,2),3);

A = [Col1, Col2];

% Percentiles for first column (peakiness)
percentiles = prctile(A(:,1), [10 20 30 40 50 60 70 80 90]); 
m = A(:,1);
vals10 = (m <= percentiles(1)); % values in the 10th percentile
vals20 = (percentiles(1) < m & m <= percentiles(2)); % values in the 20th percentile
vals30 = (percentiles(2) < m & m <= percentiles(3)); % values in the 30th percentiles
vals40 = (percentiles(3) < m & m <= percentiles(4)); % values in the 40th percentiles
vals50 = (percentiles(4) < m & m <= percentiles(5)); % values in the 50th percentiles
vals60 = (percentiles(5) < m & m <= percentiles(6)); % values in the 60th percentiles
vals70 = (percentiles(6) < m & m <= percentiles(7)); % values in the 70th percentiles
vals80 = (percentiles(7) < m & m <= percentiles(8)); % values in the 80th percentiles
vals90 = (percentiles(8) < m & m <= percentiles(9)); % values in the 90th percentiles
vals100 = (percentiles(9) < m); % values in the 100 percentiles

vals = 0.*vals10 + 0.*vals20 + 2.*vals30 + 2.*vals40 + ...
    3.*vals50 + 6.*vals60 + 12.*vals70 + 15.*vals80 + ...
    15.*vals90 + 10.*vals100;

% Store the percentile of each peakiness entry in third column of A
A = [A, vals];

% Percentiles for second column
percentiles = prctile(A(:,2), [10 20 30 40 50 60 70 80 90]); 
m = A(:,2);
vals10 = (m <= percentiles(1)); % values in the 10th percentile
vals20 = (percentiles(1) < m & m <= percentiles(2)); % values in the 20th percentile
vals30 = (percentiles(2) < m & m <= percentiles(3)); % values in the 30th percentiles
vals40 = (percentiles(3) < m & m <= percentiles(4)); % values in the 40th percentiles
vals50 = (percentiles(4) < m & m <= percentiles(5)); % values in the 50th percentiles
vals60 = (percentiles(5) < m & m <= percentiles(6)); % values in the 60th percentiles
vals70 = (percentiles(6) < m & m <= percentiles(7)); % values in the 70th percentiles
vals80 = (percentiles(7) < m & m <= percentiles(8)); % values in the 80th percentiles
vals90 = (percentiles(8) < m & m <= percentiles(9)); % values in the 90th percentiles
vals100 = (percentiles(9) < m); % values in the 100th percentiles

vals = 0.*vals10 + 2.*vals20 + 3.*vals30 + 10.*vals40 + ...
    10.*vals50 + 10.*vals60 + 5.*vals70 + 3.*vals80 + ...
    2.*vals90 + 0.*vals100;

% Store the percentile of each firing rate entry in fourth column of A
A = [A, vals];

% Make aggregate score from both columns
A = [A, A(:,end)+A(:,end-1)];

% Get the sorting index and return this (for sorting high to low)
[~,ind] = sortrows(A,size(A,2),'descend');

end
