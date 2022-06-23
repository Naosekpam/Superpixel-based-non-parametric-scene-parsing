% Summary of offine operations
RunSystemMake3D;

disp('1. Compute image descriptors for retrieval.\n');
ComputeGistDescriptors;
disp('2. Extract SIFT descriptors and their spare codings.\n');
if ~exist(['Vocabs/siftVocab' num2str(szSiftVocab) '.mat'], 'file')
    BuildSIFTVocab;
end
ComputeSIFT;
ComputeSCode;

disp('3. Compute superpixels and their ground truth labels (slow).\n');
ComputeFHSegments;
ComputeGTLabelsFH;

disp('4. Compute ground truth segments and labels for data statistics.\n');
ComputeGTSegments;
%ComputeCoocStatistics;

%disp('5. Extract superpixel descriptors.\n');
if ~exist(['Vocabs/colorVocab' num2str(szColorVocab) '.mat'])
    BuildColorVocab;
end
%ComputeSegmentDescriptors;



