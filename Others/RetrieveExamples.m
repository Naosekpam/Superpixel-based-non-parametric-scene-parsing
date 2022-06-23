% load retrieval set ##NO RE-RANKING IS DONE
%K-Refers to the number of images to be retrieved
run('E:\vlfeat-0.9.21\toolbox\vl_setup.m');%####RETRIEVAL USING GIST
GistQR = im2gist(im);

%centroid=[GistDB]
 %centroid_GIST=(GistDB);
 %kdtree_gist = vl_kdtreebuild(centroid_GIST);
 %gistQuery=(single(GistQR));
      
%[index_GIST,distance_GIST] = vl_kdtreequery(kdtree_gist,centroid_GIST,gistQuery,'NumNeighbors', 4);
%retrList = trainList(index_gist);
%gsims = GistQR'*GistDB;
%[gsims,rank] = sort(gsims, 'descend');
%rank=sort(rank);
%retrList = trainList(index_gist);
%gsims = GistQR'*GistDB;
%[gsims,rank] = sort(gsims, 'descend');
%rank=sort(rank);

%retrList1 = trainList(rank(1:K));


%retrList1 = trainList(index_GIST);

%####RETRIEVAL USING VLAD
% Load image and extract SIFT features
new_image=im;
new_image = single(rgb2gray(new_image));
[~, new_sift] = vl_sift(new_image);

% Create assignment matrix
nn = vl_kdtreequery(kdtree_vlad, centroids, single(new_sift));
assignments = zeros(64, numel(nn), 'single');
assignments(sub2ind(size(assignments), nn, 1:length(nn))) = 1;
% Encode using VLAD
new_vlad = vl_vlad(single(new_sift), centroids, assignments);
%To find out the nearest neighbour of the query image
centroid=[enc];
%centroid1=centroid;
kdtree1 = vl_kdtreebuild(centroid);
Query = new_vlad;
Query=double(Query);
q=Query;
      
[index_VLAD,distance_VLAD] = vl_kdtreequery(kdtree1,centroid,q,'NumNeighbors',K) ;
retrList2=trainList(index_VLAD);

retrList=[retrList2];
%retrList= Gist + VLAD
