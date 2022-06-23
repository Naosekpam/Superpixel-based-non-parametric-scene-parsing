%Combination of Gist, VLAD, Bag of Words and Spatial Pyramid Matching for
%scene retrieval in the Non-Parametric Scene Parsing pipeline 
%Creation date : 27th February 2019
clc;
clear all;
%Link the VLFeat Library
run('F:\Study Materials\ADBU\sem2\Matlab\ComputerVision\vlfeat-0.9.21\toolbox\vl_setup.m');
%File path where images are stored
repo = 'F:/SiftFlowDatasetOrig/Images/outdoorcategories/';
filelist = dir([repo '*.jpg']);

%###############Scene Retrieval using Gist#########################
% setting up GIST Parameters:
clear param;
param.orientationsPerScale = [8 4]; % number of orientations per scale (from HF to LF)
param.numberBlocks = 4;
param.fc_prefilt = 4;
kernelsize = 2;
Nimages=2688; %No of training images
%Pre-allocate gist:
Nfeatures = sum(param.orientationsPerScale)*param.numberBlocks^2;
gist = zeros([Nimages Nfeatures]); 

%GIST evaluation of every image in the dataset
for x=1:size(filelist, 1)
    images{x}=imread([repo,filelist(x).name]);
    gist(x, :) = LMgist(images{x}, '', param);
end
%Build a KD tree for retrieval and indexing
 centroid=[gist]
 centroids=transpose(centroid);
 kdtree = vl_kdtreebuild(centroids);
 %Query image is image with index 5
Query=imread('F:/SiftFlowDatasetOrig/Images/outdoorcategories/mountain_land188.jpg');
gistQuery(1,:) = LMgist(Query, '', param); 
q=transpose(gistQuery);
      
%Show the index and distance of the nearest neighbour n images similar to the query image
 [index_gist,distance_gist] = vl_kdtreequery(kdtree,centroids,q,'NumNeighbors', 5);
%Show the resultant retrieved images
for i=1:5
    im = imread([repo filelist(index_gist(i)).name]) ;
    figure
    imshow(im);
end


%##########Scene retrieval using VLAD #######################
%Read the images and extract SIFT descriptors from the images
for i = 1:size(filelist, 1)
    I = imread([repo filelist(i).name]) ;
    I = im2single(rgb2gray(I));
    [f,d] = vl_sift(I) ;
    sift_descr{i} = d;
end

%Build KD tree
all_descr = single([sift_descr{:}]);
centroids = vl_kmeans(all_descr, 64);
kdtree = vl_kdtreebuild(centroids);
%Initialization of the VLAD encoding. For each image VLAD encoding is
%64x128=8192D
enc = zeros(64*128, numel(sift_descr));

for k=1:numel(sift_descr)
    % Create assignment matrix
    nn = vl_kdtreequery(kdtree, centroids, single(sift_descr{k}));
    assignments = zeros(64, numel(nn), 'single');
    assignments(sub2ind(size(assignments), nn, 1:length(nn))) = 1;
    % Encode using VLAD
    enc(:, k) = vl_vlad(single(sift_descr{k}),centroids, assignments);
end

% Load image and extract SIFT features of the Query image
new_image = imread('F:/TrialDataset/coast_bea1.jpg');
new_image = single(rgb2gray(new_image));
[~, new_sift] = vl_sift(new_image);

% Create assignment matrix
nn = vl_kdtreequery(kdtree, centroids, single(new_sift));
assignments = zeros(64, numel(nn), 'single');
assignments(sub2ind(size(assignments), nn, 1:length(nn))) = 1;
% Encode using VLAD
new_vlad = vl_vlad(single(new_sift), centroids, assignments);
%To find out the nearest neighbour of the query image
centroid=[enc];
centroids=centroid;
kdtree = vl_kdtreebuild(centroids);
Query = new_vlad;
Query=double(Query);
q=Query;
%Show the index and distance of the nearest neighbour n images similar to the query image      
[index_VLAD,distance_VLAD] = vl_kdtreequery(kdtree,centroids,q,'NumNeighbors', 5) ;
for i=1:5
    im = imread([repo filelist(index_VLAD(i)).name]) ;
    figure
    imshow(im);
end

%##########Scene retrieval using Bag of Visual Words #######################
%Read the images and extract SIFT descriptors from the images
sift_descr = {};
for i = 1:size(filelist, 1)
    I = imread([repo filelist(i).name]) ;
    I = im2single(rgb2gray(I));
    [f,d] = vl_sift(I) ;
    sift_descr{i} = d;
end

%Perform clustering and build KD tree
all_descr = single([sift_descr{:}]);
centroids = vl_kmeans(all_descr, 64);
forest = vl_kdtreebuild(centroids);
%Build BoW histogram for each training images of dimension 64
for i = 1:size(filelist, 1)
  I = imread([repo filelist(i).name]) ;
  I = im2single(rgb2gray(I));
  [f,d] = vl_sift(I) ;
  [index,dist] = vl_kdtreequery(forest, centroids, single(d));
  index=double(index);
  feature_hist = hist(index, 64);
  feature_hist = feature_hist ./ sum(feature_hist);
  feature_hist = feature_hist ./ norm(feature_hist);
  image_feats(i, :) = feature_hist; 
end
%Build BoW histogram for the query image based on the codebook generated
fprintf('Now for the Query image\n');
new_image = imread('F:/TrialDataset/coast_bea1.jpg');
new_image = single(rgb2gray(new_image));
[~, new_sift] = vl_sift(new_image);

  [index_t,dist_t] = vl_kdtreequery(forest, centroids, single(new_sift));
  index_t=double(index_t);
  feature_hist_t = hist(index_t, 64);
  feature_hist_t = feature_hist_t ./ sum(feature_hist_t);
  feature_hist_t = feature_hist_t ./ norm(feature_hist_t);
  image_feats_t = feature_hist_t;
  
image_feats=transpose(image_feats);
forest_new = vl_kdtreebuild(image_feats);%kd trees based on the histogram of each training images
image_feats_t=transpose(image_feats_t);
%Show the index and distance of the nearest neighbour n images similar to the query image  
[index_BoW,dist_BoW] = vl_kdtreequery(forest_new,image_feats,image_feats_t,'NumNeighbors',5)
for i=1:5
    Im = imread([repo filelist(index_BoW(i)).name]) ;
    figure;
    imshow(Im);
end

%##########Scene retrieval using Spatial Pyramid Matching #######################
%folder for saving resultant SIFT features
data_dir = 'C:/Program Files/MATLAB/R2016a/bin/CS231BPA3-master/data/';
%Read the images and extract SIFT descriptors from the images
%filelist = dir([repo '*.jpg']);
fnames = dir(fullfile(repo, '*.jpg'));
num_files = size(fnames,1);
filenames = cell(num_files,1)

for f = 1:num_files
	filenames{f} = fnames(f).name;
end

% return pyramid descriptors for all files in filenames
pyramid_all = BuildPyramid(filenames,repo,data_dir);
pyramid_all=transpose(pyramid_all);
forest_new = vl_kdtreebuild(pyramid_all);

a=pyramid_all;
image_feats_t=a(:,5);
%Show the index and distance of the nearest neighbour n images similar to the query image  
[index_SPM,dist_SPM] = vl_kdtreequery(forest_new,pyramid_all,image_feats_t,'NumNeighbors',5);
for i=2:6
    Im = imread([repo filelist(index_SPM(i)).name]) ;
    figure;
    imshow(Im);
end

%%%########Now save all the indexes in one array without re-ranking#######
index_gist1=transpose(index_gist);
index_VLAD1=transpose(index_VLAD);
index_BoW1=transpose(index_BoW);
index_SPM1=transpose(index_SPM);

index_all=[index_gist1 index_VLAD1 index_BoW1 index_SPM1];



%############The followings are trial codes for re-ranking of the retrieval
%set##################

index_gist1=transpose(index_gist);
index_VLAD1=transpose(index_VLAD);
index=[index_gist1 index_VLAD1]

%y = zeros(size(index_gist1));
for i = 1:length(index_gist1)
gist1.index(i) = i;
gist1.content(i)=index_gist1(i);
end

for i = 1:length(index_VLAD1)
vlad1.index(i) = i;
vlad1.content(i)=index_VLAD1(i);
end

for i=1:length(index_gist1)
 for j=1:length(index_VLAD1)
   if (gist1.content(i)==vlad1.content(j))
      avg.count(i)=(gist1.index(i)+vlad1.content(j));
   else   
       avg.count(i)=0;
   end
 end   
end

%#############END OF TRIAL OF RE_RANKING#########


 