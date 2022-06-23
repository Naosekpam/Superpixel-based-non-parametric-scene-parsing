function [descrs, descrsContext] = compute_superpixel_descriptors(im, baseRegions, nBaseRegions, siftVocab, colorVocab)

% some default parameters
llcK = 5;
xLocZones = 6;
yLocZones = 6;
se = strel('diamond', 10);
%se = strel('square', 10);
%
%% sift histogram
imsize = size(im);
if length(imsize)<3
   im = cat(3,im,im,im);
end
[frames sifts] = vl_phow(single(rgb2gray(im)),'Step',4);
imsize = imsize(1:2);

[ix w] = LLCEncode(single(sifts), single(siftVocab), llcK);

%% compute superpixel descriptors
[bh, bv] = phogGradients(im);
descrs = cell(1,length(baseRegions));
bboxes = cell(1,length(baseRegions));
descrsContext = cell(1,length(baseRegions));
bboxesContext = cell(1,length(baseRegions));

 
   %%Added extra
    
for kk = 1:length(baseRegions)
      %Added extra to get the mask to get dilated RGB histogram  
       
     %[foo, borders, bb] = get_int_and_borders(mask4);
     %maskCrop = mask4(bb(1):bb(2),bb(3):bb(4));       
     
     %mask2=borders(:,:,5);
     %mask2=uint8(mask2)
        
    
        %computing sift histograms for base regions
        loc = sub2ind(imsize, round(frames(2,:)'), round(frames(1,:)'));
        regions=baseRegions{kk}(loc);
        regions = repmat(regions', [llcK 1]);
        wordIndex = (regions-1)*size(siftVocab,2)+ix;
        wordHist = accumarray(wordIndex(:),w(:),[nBaseRegions(kk)*size(siftVocab,2) 1],@sum);
        wordHist = reshape(wordHist,[size(siftVocab,2) nBaseRegions(kk)]);

        %computing RGB color histograms for base regions
        colors = reshape(shiftdim(im,2),[3 size(im,1)*size(im,2)]);
        colors = vl_ikmeanspush(colors,colorVocab);
        colorHist = accumarray([colors(:) baseRegions{kk}(:)],...
        ones(numel(colors),1), [size(colorVocab,2) nBaseRegions(kk)]);
    
     %computing HSV color  histograms for base regions ###ADDED FEATURES
    % im3=rgb2hsv(im);
    % im3=im2uint8(im3);   
     %colors_hsv = reshape(shiftdim(im3,2),[3 size(im3,1)*size(im3,2)]);
     %colors_hsv = vl_ikmeanspush(colors_hsv,colorVocab);
     %colorHist_hsv = accumarray([colors_hsv(:) baseRegions{kk}(:)],...
     %ones(numel(colors_hsv),1), [size(colorVocab,2) nBaseRegions(kk)]);
        
     %computing location histograms for base regions
        [y x] = ndgrid(1:size(im,1),1:size(im,2));
        y = ceil(y/size(im,1)*yLocZones);
        x = ceil(x/size(im,2)*xLocZones);
        locations = (y-1)*xLocZones+x;
        locHist = accumarray([locations(:) baseRegions{kk}(:)],...
        ones(numel(locations),1), [xLocZones*yLocZones nBaseRegions(kk)]);       
   
        %computing shape histograms for base regions
        bboxes{kk} = zeros(4,nBaseRegions(kk));
        for i = 1:nBaseRegions(kk)
            mask = (baseRegions{kk} == i);
            rangx = find(sum(mask,1));  
            rangy = find(sum(mask,2));
            minx = rangx(1); maxx = rangx(end);
            miny = rangy(1); maxy = rangy(end);
            bboxes{kk}(:,i) = [minx;miny;maxx;maxy];
        end
        shapeHist = phog2(bh, bv, bboxes{kk});
         
       % hogHist=extractHOGFeatures(base{kk},'CellSize',[40 40]);
        %computing dilated RGB color histograms for base regions
       % im1=im;
       % im1= im2uint8(im1);
       % im1=im1.*repmat(mask2,[1,1,3]);
        %colors1 = reshape(shiftdim(im1,2),[3 size(im1,1)*size(im1,2)]);
        %colors1 = vl_ikmeanspush(colors1,colorVocab);
        %colorHist1 = accumarray([colors1(:) baseRegions{kk}(:)],...
        %ones(numel(colors1),1), [size(colorVocab,2) nBaseRegions(kk)]);
        
        %L1-normalization
        wordHist = bsxfun(@times, wordHist, 1./(sum(wordHist)+eps));
        colorHist = bsxfun(@times,colorHist,1.0./(sum(colorHist)+1e-10));
        locHist = bsxfun(@times,locHist,1.0./(sum(locHist)+1e-10));
        shapeHist = bsxfun(@times,shapeHist,1.0./(sum(shapeHist)+1e-10));
       %colorHist1 = bsxfun(@times,colorHist1,1.0./(sum(colorHist1)+1e-10));
      % colorHist_hsv=bsxfun(@times,colorHist_hsv,1.0./(sum(colorHist_hsv)+1e-10));
        
       descrs{kk} = single([wordHist;colorHist;locHist;shapeHist]);

        
        
        %Contexual info
        wordHistContext = zeros(size(siftVocab,2), nBaseRegions(kk),'single');
        colorHistContext = zeros(size(colorVocab,2), nBaseRegions(kk),'single');
        locHistContext = zeros(xLocZones*yLocZones, nBaseRegions(kk),'single');
        bboxesContext{kk} = zeros(4,nBaseRegions(kk));
        for i = 1:nBaseRegions(kk)
            mask = (baseRegions{kk} == i);
            mask = imdilate(mask, se);
            rangx = find(sum(mask,1));  
            rangy = find(sum(mask,2));
            minx = rangx(1); maxx = rangx(end);
            miny = rangy(1); maxy = rangy(end);
            bboxesContext{kk}(:,i) = [minx;miny;maxx;maxy];
            validLoc = mask(loc);
            subIndex = ix(:,validLoc>0);
            subCoeff = w(:,validLoc>0);
            wordHistContext(:,i) = accumarray(subIndex(:),subCoeff(:),[size(siftVocab,2) 1],@max);
            subColors = colors(:,mask>0);
            colorHistContext(:,i) = accumarray(subColors(:),ones(numel(subColors),1), [size(colorVocab,2) 1]);
            subLocations = locations(mask>0);
            locHistContext(:,i) = accumarray(subLocations(:),ones(numel(subLocations),1), [xLocZones*yLocZones 1]); 
        end
        shapeHistContext = phog2(bh, bv, bboxesContext{kk});

        
        %computing dilated color histograms for base regionn   
      %  colors1_contxt = reshape(shiftdim(im1,2),[3 size(im1,1)*size(im1,2)]);
       % colors1_contxt = vl_ikmeanspush(colors1_contxt,colorVocab);
        %colorHist1_contxt = accumarray([colors1_contxt(:) baseRegions{kk}(:)],...
        %ones(numel(colors1_contxt),1), [size(colorVocab,2) nBaseRegions(kk)]);
         
       %computing color hsv histograms for base regions
       % colors_hsv_contxt = reshape(shiftdim(im3,2),[3 size(im3,1)*size(im3,2)]);
       % colors_hsv_contxt = vl_ikmeanspush(colors_hsv_contxt,colorVocab);
      %  colorHist_hsv_contxt = accumarray([colors_hsv_contxt(:) baseRegions{kk}(:)],...
      %  ones(numel(colors_hsv_contxt),1), [size(colorVocab,2) nBaseRegions(kk)]);
        %L1-normalization
        wordHistContext = bsxfun(@times, wordHistContext, 1./(sum(wordHistContext)+eps));
        colorHistContext = bsxfun(@times,colorHistContext,1.0./(sum(colorHistContext)+1e-10));
        locHistContext = bsxfun(@times,locHistContext,1.0./(sum(locHistContext)+1e-10));
        shapeHistContext = bsxfun(@times,shapeHistContext,1.0./(sum(shapeHistContext)+1e-10));
       
         %colorHist1_contxt=bsxfun(@times,colorHist1_contxt,1.0./(sum(colorHist1_contxt)+1e-10));
        % colorHist_hsv_contxt=bsxfun(@times,colorHist_hsv_contxt,1.0./(sum(colorHist_hsv_contxt)+1e-10));
         descrsContext{kk} = single([wordHistContext;colorHistContext;locHistContext;shapeHistContext]);
        
    
    
    
end