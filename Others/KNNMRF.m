
RunSystemMake3D;

visualize = true;
kernel = 'Chi2'; % kernel for comparing superpixel descriptors
%Can be change to Hellinger or Intersection kernel
K = 3; %# no of retrieval images
knn =100 ;% # of k nearest neighbors for segment classification
alpha = 6; % parameter for MRF pairwise term

%load(gistDir);
load(['Vocabs/colorVocab' num2str(szColorVocab) '.mat']);
load(['Vocabs/siftVocab' num2str(szSiftVocab) '.mat']);


for r =1:length(testList)
    
     r
    
    query = testList{r};
    
    imfile = ['../datasets/Make3D/Test/Images/' query '.jpg'];
    im = imread(imfile);
  
    RetrieveExamples1;
    %Features extraction of the retrieval set and the query image
    ComputeSegmentDescriptors1;
    
    % superpixel classification of the test super-pixels based on the the
    % NN superpixels of the Q image
    ClassifySuperpixels;

    % infer labels by MRF
    InferMRF;
    save([outputDir query '.mat'],'predictLabel');
    toc;
    
    % visualize
    if visualize
        draw_label_image(im/255, predictLabel, palette, ['unlabeled',classes]);
        print([outputDir, query, '.png'],'-dpng','-r96');
    end
        
    %disp('please enter any key to continue');

    
   %pause;                                           
    close all;
                                                       
end

