function fixations = calcRWAFixations(gazeX,gazeY,sampleRate,maxDispersion,minDuration)

% sampleRate = 200;
% maxDispersion = 10; %pixels
% minDuration = 300; %ms

% startIndex = 1;
% endIndex = 1;
% fixations = false(size(gazeX));
% for k=2:length(gazeX)
%     dispersion = max([max(abs(gazeX(startIndex:k)-gazeX(startIndex))),max(abs(gazeY(startIndex:k)-gazeY(startIndex)))]);
%     if dispersion <= maxDispersion
%         endIndex = k;
%     else
%         duration = (endIndex-startIndex+1)/sampleRate*1000;
%         if duration >= minDuration
%             fixations(startIndex:endIndex) = true;
%         end
%         startIndex = k;
%     end
% end

fixations = all(movmad(abs([gazeX,gazeY]),5) < maxDispersion,2); %calc mad over 5 samples (~25ms at 200Hz)

clust_info = bwconncomp(fixations);

clust_remove = find(cellfun(@numel,clust_info.PixelIdxList) < minDuration/1000*sampleRate);

for k=1:length(clust_remove)
    fixations(clust_info.PixelIdxList{clust_remove(k)})=0;
end