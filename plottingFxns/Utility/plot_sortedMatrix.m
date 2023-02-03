function plot_sortedMatrix(matrix,indx,direction)

A = matrix(indx,:);
C = nanmean(A(:,90:100),2) - nanmean(A(:,101:110),2);
% C = normalize(C);
C(isnan(C)) = -Inf;
[~,index] = sortrows(C,direction);
% A(isinf(A)) = NaN;
B = A(index,:);
imagesc(normalize(B))
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
end