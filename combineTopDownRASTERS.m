function combineTopDownRASTERS(session_name, size_vert, size_horiz, numrow)
% script to combine rasters into one giant png for one session
% MGC 3/1/2019
% functionalized by FKM on 7/1/19
%% params
% Input:
%   % size of images
%   size_vert = 167; % changed from 1042 to 346
%   size_horiz = 667; % changed from 333 to 667 for repeating tracks
%   numrow = 6; % number of rows in final image
% 
%   % session names
%   session_name = {'F3_190625_johnrepeatingtrack_meth2'}; % new output from singeSessionRasterplots

% where to save images
image_save_dir = '/Users/KeiMasuda/Desktop/fkm_analysis/rasters_dch/pretty_rasters/';
image_dir = fullfile('/Users/KeiMasuda/Desktop/fkm_analysis/rasters_dch/');

% image_save_dir = '/Users/KeiMasuda/Desktop/fkm_analysis/rasters/pretty_rasters/';
% image_save_dir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/E2/E2_190614_johncontrasttrack_train1_g0/E2_190614_johncontrasttrack_train1';

if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end
%% Combine images for one session
animalName = extractBefore(session_name,'_');
sessionDate = extractBefore(extractAfter(session_name,'_'),'_');

%     image_dir = fullfile('/Users/KeiMasuda/Desktop/fkm_analysis/rasters/');
% image_dir = fullfile('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/E2/E2_190614_johncontrasttrack_train1_g0/E2_190614_johncontrasttrack_train1/select_trials_pretty_rasters');

% get png file names
png_files = dir(sprintf('%s%s_%s_*.png',image_dir,animalName,sessionDate));
png_files = {png_files.name};


% create empty matrix for holding final image
numcol = ceil(numel(png_files)/numrow);
final_image = nan(numrow*size_vert,numcol*size_horiz,3);

% fill in final_image with data from png files
for i = 1:numel(png_files)
    fprintf('%s, file %d/%d\n',session_name,i,numel(png_files));

    % read image
    dat = imread(fullfile(image_dir,png_files{i}),'png');

    % get row number and column number
    row = ceil(i/numcol);
    col = mod(i-1,numcol)+1;

    % enter data into final image
    final_image((row-1)*size_vert+1:row*size_vert,(col-1)*size_horiz+1:col*size_horiz,:) = dat;
end

% convert to uint8
final_image = uint8(final_image);

% write to file
imwrite(final_image,fullfile(image_save_dir,sprintf('%s_all_rasters_combined.png',session_name)));

fprintf('Done Combining\n');
end