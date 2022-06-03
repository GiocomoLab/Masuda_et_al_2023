% PART 1 - https://www.crowdcast.io/e/sara-sollas-world-wide, 15:00-25:00
% Load synthetic data and find low-dimensional manifold.
% Project data on neural modes; see how latent trajectory generates data.
% Jacob Bakermans, February 2021
close all;
%% Section 1: Load data 
% In this section we'll generate some synthetic neural data for 3 neurons,
% and plot the data over time and as a trajectory in neuron activity space

% Load data
nCells = size(cells.spatialFR2,1);

all_fr = [];
ds_factor = 100;
trial = squeeze(cell2mat(struct2cell(cells.trial)));
trial = trial(:,1);
trial_ds = downsample(trial, ds_factor); 

% gaussian filter for smoothing
smoothSigma = 10;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);

% Calculating smoothed firing rate
for i = 1:nCells
    cellFR = cells.FRtime(i).FRtime;
    % smooth firing rate
    fr_smoothed = conv(repmat(cellFR,1,3),gauss_filter,'same');
    fr = fr_smoothed(numel(cellFR)+1:numel(cellFR)*2);
    %
    all_fr(i, :) = fr;
end 

all_fr = all_fr';
all_fr_ds = downsample(all_fr, ds_factor);
%%
X_D = all_fr_ds'; % Data matrix: N_neurons x N_timebins = 106x175359

% Plot data: firing rate for each neuron
figure(); 
% First subplot: neuron firing rate through time
subplot(1,2,1);
hold on;
for currNeuron = 1:size(X_D,1)
    plot(X_D(currNeuron, :));
end
hold off;
% Set plot layout properties
% legend('Neuron 1', 'Neuron 2', 'Neuron 3');
ylabel('Firing Rate (1/0.2 sec)');
xlabel('Time');
title('Neuron activity through time');
% Second subplot: trajectory in neural space
subplot(1,2,2);
% Plot data
plot3(X_D(1,:), X_D(2,:), X_D(3,:), 'x-');
% Set plot layout properties
xlabel('Neuron 1 activity');
ylabel('Neuron 2 activity');
zlabel('Neuron 3 activity');
view(75, 30);
grid on;
title('Trajectory in neural space');

%% Section 2: Do PCA to get neural modes/principal components 
% In this section we'll find the flat subspace/linear manifold in the data by PCA.
% This manifold is spanned by the first two principal components, or neural modes.

% A PCA is an eigendecomposition of the covariance matrix of a series of datapoints
% Remove mean firing rate for every neuron
X_D_og = X_D;
X_D = X_D - repmat(mean(X_D, 2), [1, size(X_D,2)]);
% Calculate covariance matrix between neurons
cov_dat = X_D*X_D';
% Do eigenvalue decomposition of covariance
[V, D] = eig(cov_dat);

% V has eigenvectors in columns: these are the principal components, or neural modes
% D has corresponding eigenvalues on diagonal
% The higher the eigenvalue, the more variance of activity is explained by that neural mode
% Therefore it's useful to sort modes by eigenvalue: highest eigenvalues first

% Sort eigenvectors (columns) in V by descending eigenvalue (diagonal) in D
V = sortrows([diag(D) V'],'descend');
D = diag(V(:,1));
V = V(:,2:end)';

% Plot neural trajectory and principal components
figure();
hold on;
% Plot data
plot3(X_D(1,:), X_D(2,:), X_D(3,:), 'x-');
% Plot principal components
% for currDir = 1:106
for currDir = 1:3
    quiver3(0, 0, 0, max(abs(X_D(:)))*V(1, currDir), ...
        max(abs(X_D(:)))*V(2, currDir), max(abs(X_D(:)))*V(3, currDir), 0, ...
        'LineWidth', 4);
end
% Plot plane spanned by first two principal components
fmesh(@(s,t)V(1,1)*s+V(1,2)*t, @(s,t)V(2,1)*s+V(2,2)*t, @(s,t)V(3,1)*s+V(3,2)*t, ...
    [-1, 1])
alpha(0.5);
hold off;
% Set plot layout properties
legend('Data', 'Principal component 1',...
    'Principal component 2', 'Principal component 3', 'Manifold/subpace');
xlabel('Neuron 1 activity');
ylabel('Neuron 2 activity');
zlabel('Neuron 3 activity');
xlim([-1,1]);
ylim([-1,1]);
zlim([-1,1]);
view(75, 30);
grid on;

%% Section 3: Go from neuron activity space to neural mode space
% In this section, instead of plotting the trajectory in neuron activity space, we'll move
% the trajectory to the neural mode space: the space spanned by the neural modes
% (u and v on the slides). We'll get what Solla calls latent rates. You could think about
% latent rates as generating the measured neuron activity by linear combination.

% The first two neural modes span the manifold. Select the corresponding columns from V
V_tilde = V(:, [1:6]);

% Now move trajectory in neuron activity space to neural mode space to get latent rates
L = V_tilde' * X_D;

% Plot trajectory in neural mode space
figure(); 
plot(L(1,:), L(2,:), '-x');
axis square;
xlabel('Neural mode 1 (u)');
ylabel('Neural mode 2 (v)');
title('Latent variables: data projected on manifold, in neural mode space');

% Solla's "alternative interpretation": latent rates generate neural activity
% Use linear combinations of the latent trajectory to recover original activity
X_gen = V_tilde * L;

% Plot neuron activity through time, and compare to latent activity and recovered activity
figure(); 
% First subplot: original data
subplot(3,1,1)
hold on;
for currNeuron = 1:size(X_D,1)
    plot(X_D(currNeuron, :));
end
hold off;
% Set plot layout properties
legend('Neuron 1', 'Neuron 2', 'Neuron 3');
ylabel('Normalised spike counts (1)');
xlabel('Time (bins)');
title('Neuron activity through time');
% Second subplot: latent activity
subplot(3,1,2)
hold on;
for currNeuron = 1:size(L,1)
    s = sprintf("Latent unit %d",currNeuron);
    plot(L(currNeuron, :),"DisplayName",s);
end
hold off;
% Set plot layout properties
legend();
ylabel('Normalised spike counts (1)');
xlabel('Time (bins)');
title('Latent activity through time');
% Third subplot: recovered original data from latent units
subplot(3,1,3)
hold on;
for currNeuron = 1:size(X_gen,1)
    plot(X_gen(currNeuron, :));
end
hold off;
% Set plot layout properties
legend('Neuron 1', 'Neuron 2', 'Neuron 3');
ylabel('Normalised spike counts (1)');
xlabel('Time (bins)');
title('Recovered activity through time');

%% Section 5: peek ahead: SVD
% In this section we'll get a head start with the next bit of the talk,
% where the Singular Value Decomposition (SVD) is going to show up many times.
% Turns out that's basically the same as what we've been doing so far...

% Do a singular value decomposition of the input data X_D, 
% so that X_D = U_svd * S_svd * V_svd'
[U_svd, S_svd, V_svd] = svd(X_D);

% U_svd is a matrix with the eigenvectors of X*X' in its columns;
% V_svd is a matrix with the eigenvectors of X'*X in its columns.
% That means that U_svd should be the same as our V in the previous sections,
% except for sign flips: the same eigenvectors may point in the opposite direction

% Plot them next to each other to see if that's true
figure()
% First subplot: U from SVD
subplot(1,3,1);
imagesc(U_svd);
axis square;
title('U from SVD');
% Second subplot: V from PCA
subplot(1,3,2);
imagesc(V);
axis square;
title('V from PCA');
% Third subplot: inner products between U and V. 
% All values on the diagonal should be 1 (parallel) or -1 (opposite).
subplot(1,3,3);
imagesc(U_svd'*V);
axis square;
title('Column inner products');
colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2 - https://www.crowdcast.io/e/sara-sollas-world-wide, 45:00-55:00
% Load synthetic data from two different groups of cells, 
% to simulate different recording days, then find how latent dynamics
% across days are oriented relative to each other.
% Jacob Bakermans, February 2021
close all;
%% Section 1: Prepare data
% In this section we'll generate some synthetic neural data for two different
% recording days, with three different neurons recorded on each day
% Generate synthetic data for two different days 
[days{1}, days{2}] = generate_data_2(); % Data matrix: N_neurons x N_timebins
% Load data for two different days 
days{1} = all_fr_ds'
days{2} = all_fr_ds'; % Data matrix: N_neurons x N_timebins
% Remove mean in both
for currDay = 1:length(days)
    days{currDay} = days{currDay} - ...
        repmat(mean(days{currDay}, 2), [1, size(days{currDay},2)]);
end

% Cells: assume some overlap, two cells shared across days
example_cells{1} = [1, 2, 3];
example_cells{2} = [1, 2, 4];

% Prepare data matrices for union of all neurons across days for each day
D = length(unique(cat(1, example_cells{:})));
T = size(days{1}, 2);
X_n = zeros(D, T);
X_m = zeros(D, T);
for currCell = 1:length(example_cells{1})
    X_n(example_cells{1}(currCell), :) = days{1}(currCell, :);
    X_m(example_cells{2}(currCell), :) = days{2}(currCell, :);
end

% Plot data of both days
figure();
for currDay = 1:2
    % First subplot: neuron firing rate through time
    subplot(2, 2, 1 + (currDay - 1) *2);
    hold on;
    for currNeuron = 1:length(example_cells{currDay})
        plot(days{currDay}(currNeuron, :));
    end
    hold off;
    % Set plot layout properties
    legend(['Neuron ' num2str(example_cells{currDay}(1))],...
        ['Neuron ' num2str(example_cells{currDay}(2))], ...
        ['Neuron ' num2str(example_cells{currDay}(3))]);
    ylabel('Normalised spike counts (1)');
    xlabel('Time (bins)');
    title(['Day ' num2str(currDay) ': Neuron activity through time']);
    % Second subplot: trajectory in neural space
    subplot(2, 2, 2 + (currDay - 1) *2);
    % Plot data
    plot3(days{currDay}(1, :), ...
        days{currDay}(2, :),...
        days{currDay}(3, :), 'x-');
    % Set plot layout properties
    xlabel(['Neuron ' num2str(example_cells{currDay}(1)) ' activity']);
    ylabel(['Neuron ' num2str(example_cells{currDay}(2)) ' activity']);
    zlabel(['Neuron ' num2str(example_cells{currDay}(3)) ' activity']);
    view(75, 30);
    grid on;
    title(['Trajectory in neural space of day ' num2str(currDay)]);
end

%% Section 2: singular value decomposions
% In this section, we'll use singular value decomposition to find the low-
% dimensional manifolds in the ambient dimension D for the union of neurons
% across days. We'll plot the manifolds in the original neural activity 
% spaces because they are 3D, while the full ambient space is 4D

% Do singular value decomposition for both days. U contains basis vectors in
% D dimension in its columns and provides a basis in the space of neurons.
[U_n, S_n, V_n] = svd(X_n);
[U_m, S_m, V_m] = svd(X_m);

% Now keep only a few of the basis vectors in U: the first d dimensions, 
% corresponding to leading eigenvalues (Matlab's svd function orders basis 
% vectors by eigenvalues automatically). d is the dimensionality 
% of the low-dimensional neural manifold embedded in emperical neural space.

% Set 'flat dimension' d
d = 2;

% Select first d basis vectors from U_n and U_m
U_n_tilde = U_n(:, 1:d);
U_m_tilde = U_m(:, 1:d);

% Plot matrices we've built so far
figure();
% First subplot: U_n
subplot(2,2,1);
imagesc(U_n);
axis image;
xlabel('Columns (basis vectors)');
ylabel('Rows (neural space dimensions)');
title('$U_n$', 'interpreter','latex');
% Second subplot: U_m
subplot(2,2,2);
imagesc(U_m);
axis image;
xlabel('Columns (basis vectors)');
ylabel('Rows (neural space dimensions)');
title('$U_m$', 'interpreter','latex');
% Third subplot: U_n_tilde
subplot(2,2,3);
imagesc(U_n_tilde);
axis image;
xlabel('Columns (basis vectors)');
ylabel('Rows (neural space dimensions)');
title('$\tilde{U}_n$', 'interpreter','latex');
% Fourth subplot: U_m_tilde
subplot(2,2,4);
imagesc(U_m_tilde);
axis image;
xlabel('Columns (basis vectors)');
ylabel('Rows (neural space dimensions)');
title('$\tilde{U}_m$', 'interpreter','latex');

% Plot hyperplanes in original neural space
figure();
us = {U_n_tilde, U_m_tilde};
for currDay = 1:2
    subplot(2,1,currDay);
    hold on;    
    % Plot data
    plot3(days{currDay}(1,:), days{currDay}(2,:), days{currDay}(3,:), '.-');
    % Plot basis vectors
    for currDir = 1:2
        quiver3(0, 0, 0, ...
            max(abs(days{currDay}(:)))*us{currDay}(example_cells{currDay}(1), currDir), ...
            max(abs(days{currDay}(:)))*us{currDay}(example_cells{currDay}(2), currDir), ...
            max(abs(days{currDay}(:)))*us{currDay}(example_cells{currDay}(3), currDir), 0, ...
            'LineWidth', 4);
    end
    % Plot plane spanned by first two principal directions
    fmesh(@(s,t) us{currDay}(example_cells{currDay}(1),1)*s+us{currDay}(example_cells{currDay}(1),2)*t, ...
        @(s,t) us{currDay}(example_cells{currDay}(2),1)*s+us{currDay}(example_cells{currDay}(2),2)*t, ...
        @(s,t) us{currDay}(example_cells{currDay}(3),1)*s+us{currDay}(example_cells{currDay}(3),2)*t, ...
        [-1, 1])
    alpha(0.5);
    hold off;
    % Set plot layout properties
    legend('Data', 'Basis vector 1',...
        'Basis vector 2', 'Manifold/subpace');
    xlabel(['Neuron ' num2str(example_cells{currDay}(1)) ' activity']);
    ylabel(['Neuron ' num2str(example_cells{currDay}(2)) ' activity']);
    zlabel(['Neuron ' num2str(example_cells{currDay}(3)) ' activity']);
    xlim([-1,1]);
    ylim([-1,1]);
    zlim([-1,1]);
    view(75, 30);
    grid on;
    title(['Manifold in neural space of day ' num2str(currDay)]);
end

%% Section 3: principal angles between manifolds
% In this section we'll find the angles between the manifolds of the different
% recording days using another svd, now of the inner products of basis vectors

% Calculate inner products between basis vectors in U_n_tilde and U_m_tilde
innerProducts = U_n_tilde' * U_m_tilde;

% Do another svd on the inner product between basis vectors
[U_ip, S_ip, V_ip] = svd(innerProducts);

% Calculate the angles between hyperplanes from the diagonal of S_ip, which
% holds the ordered cosines of principal angles
angles = rad2deg(acos(diag(S_ip)));
disp(['Angles between hyperplanes are '...
    num2str(angles(1), 2) ', ' num2str(angles(2), 2)]);

% Sadly we can't actually plot the hyperplanes to inspect the angles:
% they live in the D-dimensional (D=4) space of the union of neurons recorded
% on day 1 and day 2. Just to illustrate the point, make two dummy basis matrices
% containing two 3D basis vectors each, and get the SVD of the inner products.

% First basis consists of x and y vector
u = [[1, 0, 0]', [0, 1, 0]'];
% Second basis matrix is first basis, but rotated 30 degrees along x
rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
v = rotx(deg2rad(30)) * u;

% Do SVD, and get angles from ordered cosines on diagonal of S
[~, S, ~] = svd(u'*v);
uv = rad2deg(acos(diag(S)));

% Plot basis vectors and hyperplanes
figure()
hold on;
% Plot basis vectors u
for currDir = 1:2
    quiver3(0, 0, 0, ...
        u(1, currDir), ...
        u(2, currDir), ...
        u(3, currDir), 0, ...
        'LineWidth', 4);
end
% Plot basis vectors v
for currDir = 1:2
    quiver3(0, 0, 0, ...
        v(1, currDir), ...
        v(2, currDir), ...
        v(3, currDir), 0, ...
        'LineWidth', 4);
end
% Plot plane spanned by basis vectors u
fmesh(@(s,t) u(1,1)*s+u(1,2)*t, ...
    @(s,t) u(2,1)*s+u(2,2)*t, ...
    @(s,t) u(3,1)*s+u(3,2)*t, ...
    [-1, 1])
alpha(0.5);
% Plot plane spanned by basis vectors v
fmesh(@(s,t) v(1,1)*s+v(1,2)*t, ...
    @(s,t) v(2,1)*s+v(2,2)*t, ...
    @(s,t) v(3,1)*s+v(3,2)*t, ...
    [-1, 1])
alpha(0.5);
hold off;
% Set plot layout properties
view(60,30);
legend('Basis vector u1', 'Basis vector u2', 'Basis vector v1', 'Basis vector v2',...
    'Hyperplane spanned by u', 'Hyperplane spanned by v');
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
title({'Illustration of hyperplane angle between basis sets',...
    ['Angles between hyperplanes: ' num2str(uv(1)) ', ' num2str(uv(2))]});
%% PART 3 - https://www.crowdcast.io/e/sara-sollas-world-wide, 55:00-1:00:00
% Use basis vectors for different days in common ambient space to align
% latent dynamics via Canonical Correlation Analysis, using QR-decomposition
% Jacob Bakermans, February 2021
close all;
%% Section 1: Prepare data (copy from part 2)
% In this section, we'll repeat the steps from part 2: load data, combine
% into space of union of neurons across days, do svd, and find basis for 
% low-dimensional manifold for each day

% Generate synthetic data for two different days 
[days{1}, days{2}] = generate_data_2(); % Data matrix: N_neurons x N_timebins
% Remove mean in both
for currDay = 1:length(days)
    days{currDay} = days{currDay} - ...
        repmat(mean(days{currDay}, 2), [1, size(days{currDay},2)]);
end

% Cells: assume some overlap, two cells shared across days
example_cells{1} = [1, 2, 3];
example_cells{2} = [1, 2, 4];

% Prepare data matrices for union of all neurons across days for each day
D = length(unique(cat(1, example_cells{:})));
T = size(days{1}, 2);
X_n = zeros(D, T);
X_m = zeros(D, T);
for currCell = 1:length(example_cells{1})
    X_n(example_cells{1}(currCell), :) = days{1}(currCell, :);
    X_m(example_cells{2}(currCell), :) = days{2}(currCell, :);
end

% Do singular value decomposition for both days. U contains basis vectors in
% D dimensions in its columns and provides a basis in the space of neurons.
[U_n, S_n, V_n] = svd(X_n);
[U_m, S_m, V_m] = svd(X_m);

% Set 'flat dimension' d
d = 2;

% Select first d basis vectors from U_n and U_m
U_n_tilde = U_n(:, 1:d);
U_m_tilde = U_m(:, 1:d);

%% Section 2: project data on manifolds to get latent variables
% In this section, we'll project the data of both days on the low-dimensional
% manifold for the corresponding day, to get trajectories in space of latent
% variables/neural modes. These are not yet aligned.

% Project data on d-dimensional manifolds to get latent variables
L_n = U_n_tilde' * X_n;
L_m = U_m_tilde' * X_m;

% Calculate unaligned canonical correlations from inner products
CC_unaligned = L_n*L_m';

% Plot resulting latent trajectories
% Plot hyperplanes in original neural space
figure();
us = {U_n_tilde, U_m_tilde};
ls = {L_n, L_m};
for currDay = 1:2
    % First subplot for each day: trajectory in original neural activity space
    subplot(2,2,(currDay-1)*2+1);    
    hold on;    
    % Plot data
    plot3(days{currDay}(1,:), days{currDay}(2,:), days{currDay}(3,:), 'x-');
    % Plot basis vectors
    for currDir = 1:2
        quiver3(0, 0, 0, ...
            max(abs(days{currDay}(:)))*us{currDay}(example_cells{currDay}(1), currDir), ...
            max(abs(days{currDay}(:)))*us{currDay}(example_cells{currDay}(2), currDir), ...
            max(abs(days{currDay}(:)))*us{currDay}(example_cells{currDay}(3), currDir), 0, ...
            'LineWidth', 4);
    end
    % Plot plane spanned by first two principal directions
    fmesh(@(s,t) us{currDay}(example_cells{currDay}(1),1)*s+us{currDay}(example_cells{currDay}(1),2)*t, ...
        @(s,t) us{currDay}(example_cells{currDay}(2),1)*s+us{currDay}(example_cells{currDay}(2),2)*t, ...
        @(s,t) us{currDay}(example_cells{currDay}(3),1)*s+us{currDay}(example_cells{currDay}(3),2)*t, ...
        [-1, 1])
    alpha(0.5);
    hold off;
    % Set plot layout properties
    legend('Data', 'Basis vector 1',...
        'Basis vector 2', 'Manifold/subpace');
    xlabel(['Neuron ' num2str(example_cells{currDay}(1)) ' activity']);
    ylabel(['Neuron ' num2str(example_cells{currDay}(2)) ' activity']);
    zlabel(['Neuron ' num2str(example_cells{currDay}(3)) ' activity']);
    xlim([-1,1]);
    ylim([-1,1]);
    zlim([-1,1]);
    view(75, 30);
    grid on;
    title(['Manifold in neural space of day ' num2str(currDay)]);
    % Second subplot: latent variables/trajectory in neural modes space
    subplot(2,2,(currDay-1)*2+2);
    plot(ls{currDay}(1,:), ls{currDay}(2,:), '-x');
    % Set plot layout properties
    xlabel('Neural mode 1');
    ylabel('Neural mode 2');
    title(['Latent variables on day ' num2str(currDay)]);
end

%% Section 3: Canonical Correlation Analysis (CCA)
% In this section we'll do canonical correlation analysis by QR decomposition
% of the latent variables, and svd of the orthonormal latent variables.

% If you transpose the latent variables, you'll get a T x d data matrix, 
% with the latent activity at every tim bin in the rows. 
% The QR decomposition of L' finds an orthonormal basis for L', 
% where Q holds d T-dimensional orthonormal vectors and R the coefficients
% of the linear combination that will recover input L'.

% Do the QR decomposition of transposed latent variables. Hint: to get T x d matrix 
% of orthornormal basis vectors, do "economy size" QR decomposition
[Q_n, R_n] = qr(L_n', 0);
[Q_m, R_m] = qr(L_m', 0);

% Do SVD of inner product of Qs to get aligned canonical correlation coefficients
[U_q, S_q, V_q] = svd(Q_n' * Q_m);

% Q being orthornomal means that the inner product of columns i and j of Q 
% is 1 if i == j, and 0 if i ~= j. Also L' should be recoverd by Q*R.
% Make a plot to verify both are true.
figure()
qs = {Q_n, Q_m};
rs = {R_n, R_m};
for currDay = 1:2
    % First subplot for each day: latent variables  
    subplot(2,3,(currDay-1)*3+1);    
    plot(ls{currDay}(1,:), ls{currDay}(2,:), '-x');
    % Set plot layout properties
    xlabel('Neural mode 1');
    ylabel('Neural mode 2');
    title(['Latent variables on day ' num2str(currDay)]);
    % Second subplot for each day: QR decomposition of latent variables
    % Expect inner products to be identity matrix, and different to be 0
    subplot(2,3,(currDay-1)*3+2);    
    imagesc(Q_n' * Q_n);
    % Set plot layout properties    
    xlabel('Basis vectors');
    ylabel('Basis vectors');
    title({'Inner product between columns of Q',...
        ['Maximum of L_n^T - Q_n*R_n: ' ...
        num2str(max(max(abs((Q_n * R_n) - L_n'))),2)]});
    % Second subplot for each day: latent variables in orthonormal basis 
    subplot(2,3,(currDay-1)*3+3);    
    plot(qs{currDay}(:,1), qs{currDay}(:,2), '-x');
    % Set plot layout properties
    xlabel('Orthonormal mode 1');
    ylabel('Orthonormal mode 2');
    title(['Orthonormal latent variables on day ' num2str(currDay)]);  
end

%% Section 4: Align latent variables
% In this section we'll do the final step in aligning latent dynamics: we'll
% use the outcomes of CCA to calculate and apply the linear transformation
% that takes the latent variables in the neural modes of different days into
% a single aligned manifold.

% Calculate transformation matrices M that move latent variables to manifold 
% directions in which latent variables are maximally correlated across days
M_n = inv(R_n) * U_q;
M_m = inv(R_m) * V_q;

% Use transformation matrices M to carry out the linear transformation of 
% latent variables into the aligned space
L_n_aligned = M_n' * L_n;
L_m_aligned = M_m' * L_m;

% Finally: plot results of aligned latent variables in space of maximally
% correlated neural modes found via canonical correlation analysis
figure()
for currDay = 1:2
    % First subplot for each day: latent variables  
    subplot(2,4,(currDay-1)*4+1);    
    plot(ls{currDay}(1,:), ls{currDay}(2,:), '-x');
    % Set plot layout properties
    xlabel('Neural mode 1');
    ylabel('Neural mode 2');
    title(['Latent variables on day ' num2str(currDay)]);
    % Second subplot for each day: latent variables in orthonormal basis 
    subplot(2,4,(currDay-1)*4+2);    
    plot(qs{currDay}(:,1), qs{currDay}(:,2), '-x');
    % Set plot layout properties
    xlabel('Orthonormal mode 1');
    ylabel('Orthonormal mode 2');
    title(['Orthonormal latent variables on day ' num2str(currDay)]);  
end
% One big combined subplot: both latent variables in the same aligned space
subplot(2,4,[3,4,7,8]);
hold on;
plot(L_n_aligned(1,:), L_n_aligned(2,:), '-x');
plot(L_m_aligned(1,:), L_m_aligned(2,:), '-x');
hold off;
% Set plot layout properties
legend('Day 1', 'Day 2');
xlabel('Aligned mode 1');
ylabel('Aligned mode 2');
title('Aligned latent variables across days');

%% Section 5: Post scriptum - was QR decomposition necessary?
% In this section we'll arrive at the same result taking an alternative
% route: directly use the leading factors from the initial SVD of the data,
% instead of calculating the QR decomposition of the latent variables. The 
% factors (columns of V in the initial SVD; principal components across time)
% are orthonormal already. What is the advantage of doing QR decomposition?

% Select first d basis vectors from V_n and V_m. These are 'factors', 
% vectors of linear combinations across time bins instead of neurons.
V_n_tilde = V_n(:, 1:d);
V_m_tilde = V_m(:, 1:d);

% Do SVD of inner products between factors
[U_v, S_v, V_v] = svd(V_n_tilde' * V_m_tilde);

% Now the matrix of singular values of the data is going to play to role of
% R in the QR decomposition: it holds coefficients for linear combination of
% basis vectors. Only select the first d rows and columns.
S_n_tilde = S_n(1:d, 1:d);
S_m_tilde = S_m(1:d, 1:d);

% Calculate matrices that transform data to aligned neural mode space -
% equivalent to the matrices M_m and M_n before - from SVD results
M_n_svd = inv(S_n_tilde) * U_v;
M_m_svd = inv(S_m_tilde) * V_v;

% Apply transformation to get latent variables in aligned space
L_n_aligned_svd = M_n_svd' * L_n;
L_m_aligned_svd = M_m_svd' * L_m;

% Make a plot to compare the transformation matrices from QR decomposition
% to those obtained directly from factors of original SVD
figure();
% First subplot: M_n
subplot(2,2,1);
imagesc(M_n);
axis image;
xlabel('Columns');
ylabel('Rows');
title('M_n from QR decomposition');
% Second subplot: M_m
subplot(2,2,2);
imagesc(M_m);
axis image;
xlabel('Columns');
ylabel('Rows');
title('M_m from QR decomposition');
% Third subplot: M_n_svd
subplot(2,2,3);
imagesc(M_n_svd);
axis image;
xlabel('Columns');
ylabel('Rows');
title('M_n directly from SVD');
% Fourth subplot: U_m_tilde
subplot(2,2,4);
imagesc(M_m_svd);
axis image;
xlabel('Columns');
ylabel('Rows');
title('M_m directly from SVD');

% Now plot basis vectors for alignment and aligned latent variables for
% both days, on different rows for different approaches, to compare them
figure()
vs = {V_n_tilde, V_m_tilde};
for currDay = 1:2
    % First row subplots: orthonormal latent variables for QR decomposition
    subplot(2,3,currDay);    
    plot(qs{currDay}(:,1), qs{currDay}(:,2), '-x');
    % Set plot layout properties
    xlabel('Orthonormal mode 1');
    ylabel('Orthonormal mode 2');
    title(['Orthonormal latent variables Q on day ' num2str(currDay)]); 
    % Second row subplots: factors directly from initial SVD
    subplot(2,3,3 + currDay);    
    plot(vs{currDay}(:,1), vs{currDay}(:,2), '-x');
    % Set plot layout properties
    xlabel('Factor 1');
    ylabel('Factor 2');
    title(['Factors (PCs across time) of V on day ' num2str(currDay)]);  
end
% Fourth subplot for each day: aligned latent variables using QR decomposition
subplot(2,3,3);
hold on;
plot(L_n_aligned(1,:), L_n_aligned(2,:), '-x');
plot(L_m_aligned(1,:), L_m_aligned(2,:), '-x');
hold off;
% Set plot layout properties
legend('Day 1', 'Day 2');
xlabel('Aligned mode 1');
ylabel('Aligned mode 2');
title('Aligned latent from QR decomposition');
% Fourth subplot for each day: aligned latent variables using initial SVD factors
subplot(2,3,6);
hold on;
plot(L_n_aligned_svd(1,:), L_n_aligned_svd(2,:), '-x');
plot(L_m_aligned_svd(1,:), L_m_aligned_svd(2,:), '-x');
hold off;
% Set plot layout properties
legend('Day 1', 'Day 2');
xlabel('Aligned mode 1');
ylabel('Aligned mode 2');
title('Aligned latent from initial SVD factors');
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
title({'Illustration of hyperplane angle between basis sets',...
    ['Angles between hyperplanes: ' num2str(uv(1)) ', ' num2str(uv(2))]});

