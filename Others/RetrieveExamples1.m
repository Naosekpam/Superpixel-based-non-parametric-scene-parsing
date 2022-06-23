% Computer retrieval set of a given query without re-ranking using Gist,
% VLAD and SPM
%K-Refers to the number of images to be retrieved, so if K=3, total number
%numberof retrieval set is 3x3=9
%of retriev
imDir_test = 'E:\SceneParsibg\CorrectPrediction\GTObtain\testList\test\';
run('E:\vlfeat-0.9.21\toolbox\vl_setup.m');
%-------------RETRIEVAL USING GIST-------------------
%GistQR = im2gist(im);
%gsims = GistQR'*GistDB;
%[gsims,rank] = sort(gsims, 'descend');
%retrList_Gist = trainList(rank(1:2));

%-------------RETRIEVAL USING VLAD---------------------
% Load image and extract SIFT features
new_image=im;
new_image = single(rgb2gray(new_image));
[~, new_sift] = vl_sift(new_image);
% Create assignment matrix
nn = vl_kdtreequery(kdtree_vlad, centroids, single(new_sift));
assignments = zeros(128, numel(nn), 'single');
assignments(sub2ind(size(assignments), nn, 1:length(nn))) = 1;
% Encode using VLAD
new_vlad = vl_vlad(single(new_sift), centroids, assignments);
%To find out the nearest neighbour of the query image
centroid=[enc];
kdtree1 = vl_kdtreebuild(centroid);
Query = new_vlad;
Query=double(Query);
q=Query;
[index_VLAD,distance_VLAD] = vl_kdtreequery(kdtree1,centroid,q,'NumNeighbors',1) ;
retrList_VLAD=trainList(index_VLAD);
%-------------RETRIEVAL USING SPM---------------------

%pyramid_test1=pyramid_test(r,:);
%pyramid_test1=transpose(pyramid_test1);

%pyramid_test1=transpose(pyramid_test1);
%[index_SPM,dist_SPM] = vl_kdtreequery(forest_SPM,pyramid_train,pyramid_test1,'NumNeighbors',2);
%retrList_SPM=trainList(index_SPM);

retrList=[retrList_VLAD];


