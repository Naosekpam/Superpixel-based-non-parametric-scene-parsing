clear all;
close all;
RunSystemMake3D;
stage = 1;

for i = 1:length(testList)
    
    i
    query = testList{i};
    imfile = ['C:/Program Files/MATLAB/R2016a/bin/Yang/datasets/Make3D/Test/Images/' query '.jpg'];
    im = imread(imfile);
    
    % convert labels
    load(['C:/Program Files/MATLAB/R2016a/bin/Yang/output/Make3D/results/' query  '.mat']);
    draw_label_image(im/255,predictLabel, palette, ['unlabeled',classes]);
    labelstr = ['C:/Program Files/MATLAB/R2016a/bin/Yang/output/Make3D/results/' query   '.png'];
    print(labelstr,'-dpng','-r96');
    
         figure; imshow(im);     

     pause;
    close all;
    
end

