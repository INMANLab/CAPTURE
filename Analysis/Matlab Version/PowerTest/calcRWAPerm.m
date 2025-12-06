function PM = calcRWAPerm(pwr,bpwr,trialtype)
%pwr (time x freq x trial), bpwr (1 x freq) -> baseline power, trialtype is
%-1 or 1 for calculating difference between specgrams or 0 for a single
%specgram.

ntime = size(pwr,1);
nfreq = size(pwr,2);
ntrials = size(pwr,3);
ntrials1 = sum(trialtype==-1);
ntrials2 = sum(trialtype==1);
nperm = 1000;

if isscalar(trialtype)
    permtype = 1; %shift
else
    permtype = 0; %ttest
end

rng('shuffle'); %seed the random stream with clock time
PM = nan(ntime,nfreq,nperm);
% wait_msg = parfor_wait(nperm);
parfor (m=1:nperm,8) %permutations
%     wait_msg.Send;
    if permtype %shift
        cutpoint = randi(ntime,[1,ntrials]);
        pwr_shift = zeros(size(pwr));
        for k=1:ntrials
            pwr_shift(:,:,k) = circshift(pwr(:,:,k),cutpoint(k),1); %time x freq x trial (shift time)
        end
        PM(:,:,m) = 10*log10(mean(pwr_shift,3,"omitnan")./bpwr); %time x freq;
    else %ttest
        ttshuf = trialtype(randperm(ntrials));
        tnum = squeeze(mean(pwr(:,:,ttshuf==-1),3,"omitnan")-mean(pwr(:,:,ttshuf==1),3,"omitnan")); %trial avg diff (time x freq)
        tdenom = sqrt(std(pwr(:,:,ttshuf==-1),0,3,"omitnan").^2./ntrials1+std(pwr(:,:,ttshuf==1),0,3,"omitnan").^2./ntrials2);
        PM(:,:,m) = tnum./tdenom; %time x freq;
    end
end
% wait_msg.Destroy;