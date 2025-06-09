function combineGoProVids(dir_gopro, varargin)

% dir_gopro = "\\rolstonserver\D\Data\Real World Navigation Cory\RW2\Original\Walk1\Gopro";

d = dir(dir_gopro);
d(~contains({d.name},'mp4','ignorecase',true)) = [];

pt = regexp(dir_gopro,'RW\d','match','once');
wk = regexp(dir_gopro,'Walk\d','match','once');

fname_tmp = fullfile(dir_gopro,'tmp.txt');
fname_comb = fullfile(dir_gopro,sprintf('Gopro_%s_%s_Combined.MP4',pt,wk));

fid = fopen(fname_tmp,'w');
for k=1:length(d)
    str = fullfile(d(k).folder,d(k).name);
    str = regexprep(str,'\\\','\\\\');
    str = regexprep(str,'\\','\\\');
    str = regexprep(str,'\s','\\ ');
    fprintf(fid,'file %s\r\n',str);
end
fclose(fid);

cmd = sprintf('ffmpeg -f concat -safe 0 -i "%s" -c copy "%s"',fname_tmp,fname_comb);

[status,cmdout] = system(cmd);

delete(fname_tmp);

if status
    disp('Failure.')
else
    disp('Success!')
end

