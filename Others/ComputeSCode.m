% Compute sparse coding using LLC
N = length(trainList);

llcK = 5;
load(['Vocabs/siftVocab' num2str(szSiftVocab) '.mat']);

tic
for f = 1:N

    name_segment = trainList{f};
    imfile = [imDir name_segment '.jpg'];
    outputfile = [scodeDir name_segment '.mat'];
    if exist(outputfile,'file')
        continue;
    end
    load([siftDir name_segment '.mat']);
    [ix w] = LLCEncode(single(sifts), single(siftVocab), llcK);
    
    save(outputfile, 'ix', 'w');
    if mod(f,100)==0
        fprintf('LLC: %d image in %f seconds.\n', f, toc);
    end
end
