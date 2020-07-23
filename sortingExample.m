A = rand(20,4)

percentiles = prctile(A(:,1), [10 20 30 40 50 60 70 80 90]); 
m = A(:,1);
vals10 = (m <= percentiles(1)); % values in the 10th percentile
vals20 = (percentiles(1) < m & m <= percentiles(2)); % values in the 20th percentile
vals30 = (percentiles(2) < m & m <= percentiles(3)); % values in the 30th percentiles
vals40 = (percentiles(3) < m & m <= percentiles(4)); % values in the 30th percentiles
vals50 = (percentiles(4) < m & m <= percentiles(5)); % values in the 30th percentiles
vals60 = (percentiles(5) < m & m <= percentiles(6)); % values in the 30th percentiles
vals70 = (percentiles(6) < m & m <= percentiles(7)); % values in the 30th percentiles
vals80 = (percentiles(7) < m & m <= percentiles(8)); % values in the 30th percentiles
vals90 = (percentiles(8) < m & m <= percentiles(9)); % values in the 30th percentiles
vals100 = (percentiles(9) < m); % values in the 30th percentiles

vals = vals10 + 2.*vals20 + 3.*vals30 + 4.*vals40 + ...
    5.*vals50 + 6.*vals60 + 7.*vals70 + 8.*vals80 + ...
    9.*vals90 + 10.*vals100;

A = [A, vals];

percentiles = prctile(A(:,2), [10 20 30 40 50 60 70 80 90]); 
m = A(:,2);
vals10 = (m <= percentiles(1)); % values in the 10th percentile
vals20 = (percentiles(1) < m & m <= percentiles(2)); % values in the 20th percentile
vals30 = (percentiles(2) < m & m <= percentiles(3)); % values in the 30th percentiles
vals40 = (percentiles(3) < m & m <= percentiles(4)); % values in the 30th percentiles
vals50 = (percentiles(4) < m & m <= percentiles(5)); % values in the 30th percentiles
vals60 = (percentiles(5) < m & m <= percentiles(6)); % values in the 30th percentiles
vals70 = (percentiles(6) < m & m <= percentiles(7)); % values in the 30th percentiles
vals80 = (percentiles(7) < m & m <= percentiles(8)); % values in the 30th percentiles
vals90 = (percentiles(8) < m & m <= percentiles(9)); % values in the 30th percentiles
vals100 = (percentiles(9) < m); % values in the 30th percentiles

vals = vals10 + 2.*vals20 + 3.*vals30 + 4.*vals40 + ...
    5.*vals50 + 6.*vals60 + 7.*vals70 + 8.*vals80 + ...
    9.*vals90 + 10.*vals100;

A = [A, vals];

A = [A, A(:,end)+A(:,end-1)];

[Asort,ind] = sortrows(A,size(A,2),'descend')
