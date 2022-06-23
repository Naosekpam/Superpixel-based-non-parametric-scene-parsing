%Author : Veronica Naosekpam
%Date : 5th April 2019
%First program to be run to calculate the VLAD of the training images
close all;
clear all;
%Path for the vl_setup file
run('E:\vlfeat-0.9.21\toolbox\vl_setup.m');
%Path to the training images folder
%Please change according to the system
repo = 'E:\SceneParsibg\CorrectPrediction\GTObtain\testList\TrainList\';
filelist = dir([repo '*.jpg']);
%sift_descr = {};



for i = 1:size(filelist, 1)
    I = imread([repo filelist(i).name]) ;
    I = im2single(rgb2gray(I));
    [f,d] = vl_sift(I) ;%Extract SIFT features from the training images
    sift_descr{i} = d;
end


all_descr = single([sift_descr{:}]);
centroids = vl_kmeans(all_descr, 128);
kdtree_vlad = vl_kdtreebuild(centroids);
enc = zeros(128*128, numel(sift_descr));


for k=1:numel(sift_descr)
    % Create assignment matrix
    nn = vl_kdtreequery(kdtree_vlad, centroids, single(sift_descr{k}));
    assignments = zeros(128, numel(nn), 'single');
    assignments(sub2ind(size(assignments), nn, 1:length(nn))) = 1;
    % Encode using VLAD
    enc(:, k) = vl_vlad(single(sift_descr{k}),centroids, assignments);
end

%enc stores the VLAD encoding of the trainin images.
