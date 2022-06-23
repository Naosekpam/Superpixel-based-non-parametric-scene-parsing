%Author : Veronica Naosekpam
%Date : 23rd April 2019
%Second program to be run to calculate the VLAD of the training images

%Path for the vl_setup file
run('E:\vlfeat-0.9.21\toolbox\vl_setup.m');
disp('2. Compute Spatial Pyramid of quantized SIFT descriptors for retrieval.\n');

%Path to the training/test images folder
%Please change according to the system
imDir_train = 'C:/Program Files/MATLAB/R2019a/bin/SceneParsibg/datasets/Make3D/Training/Images/';
imDir_test = 'C:/Program Files/MATLAB/R2019a/bin/SceneParsibg/datasets/Make3D/Test/Images/';

trainList_SPM =  textread('../datasets/Make3D/trainList.txt','%s');
for f = 1:length(trainList_SPM)
    trainList_SPM{f} = trainList_SPM{f};
end
%Folder path to save the SIFT descriptors computed grid-wise 
data_dir = 'C:/Program Files/MATLAB/R2019a/bin/SceneParsibg/code/data/';
%Computation of SP of the training set and build a kd tree based on
%training set
pyramid_train = BuildPyramid(trainList_SPM,imDir_train,data_dir);
%pyramid_train=transpose(pyramid_train);
pyramid_train=transpose(pyramid_train);
forest_SPM = vl_kdtreebuild(pyramid_train);
testList_SPM =  textread('../datasets/Make3D/testList.txt','%s');
pyramid_test = BuildPyramid(testList_SPM,imDir_test,data_dir);

