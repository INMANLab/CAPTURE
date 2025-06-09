function [ntp,D,Fs] = findBiopacNTP(Files,varargin)
%Data from "Biopac" folder. 

[ntp_accel,D_accel] = findChestAccelNTP(Files);

fid = fopen(Files.biopac_header_file);
hdr = fread(fid,[1,Inf],'*char');
fclose(fid);

hdr = hdr(regexp(hdr,'X,'):end-1);
hdr = regexp(regexprep(hdr,'\s+',''),',','split');
if contains(hdr(end),'}')
    a = regexp(hdr(end),'}','split','once');
    hdr(end) = a{:}(1);
end

biopac_tbl = readtable(Files.biopac_data_file);
biopac_tbl.Properties.VariableNames = hdr;

pt = str2double(regexp(regexp(Files.walk_dir,'RW\d','match','once'),'\d$','match','once'));
wk = str2double(regexp(regexp(Files.walk_dir,'Walk\d','match','once'),'\d$','match','once'));

% figure;
% subplot(2,1,1)
% plot(vecnorm(D_accel,2,2))
% title([pt,wk])
% subplot(2,1,2)
% plot(vecnorm([biopac_tbl.X,biopac_tbl.Y,biopac_tbl.Z],2,2))

Sync = cell(5,8); %pt x wk
Sync{1,1} = [1167,28515]; %chest accel x biopac accel (samples)
Sync{1,2} = [698,37344];
Sync{1,3} = [3269,44106];
Sync{1,4} = [181,20013];
Sync{1,5} = []; %no chest phone data
Sync{1,6} = [535,25695];
Sync{1,7} = [1020,31303];
Sync{2,1} = [639,14262];
Sync{2,2} = [349,11266];
Sync{2,3} = [431,16167];
Sync{2,4} = [155,6457];
Sync{2,5} = [168,7175];
Sync{2,6} = [765,14423];
Sync{2,7} = [16885,179624];
Sync{3,1} = [537,31893];
Sync{3,2} = [20177,215945];
Sync{3,3} = [602,25989];
Sync{3,4} = [318,20583];
Sync{3,5} = [248,25390];
Sync{3,6} = [544,8464];
Sync{3,7} = [269,17128];
Sync{3,8} = []; %no chest phone data
Sync{4,1} = [95,14162];
Sync{4,2} = [194,8415];
Sync{4,3} = [521,11834];
Sync{4,4} = [306,10204];
Sync{4,5} = [378,9521];
Sync{4,6} = [238,41897];
Sync{4,7} = [459,9387];
Sync{4,8} = [835,14542];
Sync{5,1} = [1249,3669];
Sync{5,2} = []; %no chest drift table
Sync{5,3} = [1226,7888];
Sync{5,4} = [829,3244];
Sync{5,5} = [1564,8256];
Sync{5,6} = []; %no chest drift table
Sync{5,7} = [1370,3665];
Sync{5,8} = []; %no chest drift or accel data


FsFull = 1000; %fixed for biopac
Fs = 100; %downsample for now otherwise eats up disk space

ntp = ((1:size(biopac_tbl,1))'-Sync{pt,wk}(2))./FsFull + ntp_accel(Sync{pt,wk}(1));
ntp = ntp(1:FsFull/Fs:end);

bNames = biopac_tbl.Properties.VariableNames;
D = nan(length(ntp),size(biopac_tbl,2));
for k=1:length(bNames)
    D(:,k) = resample(biopac_tbl.(bNames{k}),Fs,FsFull);
end
D = array2table(D,'VariableNames',bNames);



