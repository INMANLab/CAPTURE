function checkFooofModelFit(MultTrans,boxcomprng)
%This function can be added to plotMultTransSpecGramPerm to inpect the
%fooof model fits of the selected data. It uses the same parameters as the
%full walk fooof calculation

MT = MultTrans.MT;
MTd = MultTrans.MTd;

d = MTd.d_np; 
nan_idx = isnan(d); d(nan_idx) = 0;
tsec = MTd.tsec_np;
winsegsec = 2;
fs = 250;

ntime = size(d,1);
winsamp = winsegsec*fs; %window size in samples (2sec@250Hz -> 500)
stepsamp = winsamp/10; %number of samples to step (50)
win = -winsamp/2:winsamp/2; %zero-centered vector of window indices
stp = winsamp/2+1:stepsamp:ntime; %step indices adjusted so 1st sample is 1 when added to win
T = win+stp'; %matrix of indices for each window through time (windows x time indices)
rm_idx = any(T>ntime,2); %removing any indices that exceed length of data
stp(rm_idx) = [];
T(rm_idx,:) = [];
nseg = size(T,1); %number of windows
ntrials = size(d,2);

%multitaper params
cparams.Fs = fs; %250
cparams.trialave = 0;
cparams.tapers = [5,9];
cparams.fpass = [2,120];
cparams.pad = 0;
[~,f] = mtspectrumc(d(T(1,:),:),cparams);
freqs = f(f<=85); %same as fr.freqs in fooof output
fcnt = length(freqs);

%fooof params
fparams.max_n_peaks = 3; %aperiodic component has a better fit if more peaks are removed
% fparams.aperiodic_mode = 'fixed'; %'knee' (testing this 20250721), 'fixed' (default)
fparams.peak_width_limits = [1,12];
flabels = {'peak_freq','peak_height','peak_width','aperiodic_offset','aperiodic_exponent','fit_rsquared','fit_error'};

FR = nan(nseg,length(flabels),ntrials,2); %num windows x params x chan
FR_powerspec = nan(nseg,fcnt,ntrials,2); %capture power_spectrum/fooofed_spectrum/ap_fit
FR_fooofedspec = nan(nseg,fcnt,ntrials,2);
FR_apfit = nan(nseg,fcnt,ntrials,2);
for k=1:2 %fixed, knee
    if k==1
        fparams.aperiodic_mode = 'fixed';
    else
        fparams.aperiodic_mode = 'knee';
    end
    parfor m=1:nseg %by segment
        [S,f] = mtspectrumc(d(T(m,:),:),cparams); S = S'; %chan x freq
        for n=1:ntrials %by trial
            if ~any(nan_idx(T(m,:),n))
                fr = fooof(f,S(n,:),[2,85],fparams,1); %fr.peak_params = [freq, height (aperiodic removed), width] (85Hz matches upper lim in plotMultTransSpecGramPerm)
                if isempty(fr.peak_params)
                    pkp = nan(1,3);
                else
                    [~,midx] = max(fr.peak_params(:,2)); %find the highest peak
                    pkp = fr.peak_params(midx,:);
                end
                switch fparams.aperiodic_mode
                    case 'fixed'
                        FR(m,:,n,k) = [pkp,fr.aperiodic_params,fr.r_squared,fr.error];
                    case 'knee'
                        FR(m,:,n,k) = [pkp,fr.aperiodic_params([1,3]),fr.r_squared,fr.error];
                end
                FR_powerspec(m,:,n,k) = fr.power_spectrum;
                FR_fooofedspec(m,:,n,k) = fr.fooofed_spectrum;
                FR_apfit(m,:,n,k) = fr.ap_fit;
            end
        end
    end
end

plottype = 1; %1=fixed, 2=knee

%plotting
[uPtCh,~,uPtChIdx] = unique(MT(:,{'patient','regionnum'}));
mrk_style = {'o','+','*','x'}; %by region -> 1=AntHipp, 2=LatTemp, 3=Ent+Peri, 4=PostHipp+Para
mrk_size = [6,10,10,10];
mrk_color = {'r','g','b','c','m'}; %by patient -> 1:5
ptnum = uPtCh.patient;
rgnum = uPtCh.regionnum;

d_ex = squeeze(FR(:,5,:,plottype)); %grabbing the exponent
tsec_ex = tsec(stp);

tidx_box1_ex = tsec_ex>boxcomprng(1,1) & tsec_ex<boxcomprng(1,2);
tidx_box2_ex = tsec_ex>boxcomprng(2,1) & tsec_ex<boxcomprng(2,2);

d_box1_ex = median(d_ex(tidx_box1_ex,:),1);
d_box2_ex = median(d_ex(tidx_box2_ex,:),1);

nFP_box12 = nan(ntrials,2); %mean normalized across box1/box2 by patient/chan
mFP_box12 = nan(size(uPtCh,1),2);
for k=1:size(uPtCh,1)
    idx = (uPtChIdx==k);
    fp_box12 = [d_box1_ex(idx);d_box2_ex(idx)]';
    nfp_box12 = fp_box12./mean(fp_box12,2);
    nFP_box12(idx,:) = nfp_box12;
    mFP_box12(k,:) = mean(nfp_box12,1);
end


fH = figure('Position',[50,50,1850,550]);
aH = [];
aH(1) = subplot(1,3,1,'parent',fH);
boxplot(aH(1),nFP_box12,{'box1 (-.)','box2 (--)'})
[h,p,ci,stats] = ttest2(nFP_box12(:,1),nFP_box12(:,2));
ylabel(aH(1),'Normalized Aperiodic Exponent (Fooof)')
title(aH(1),sprintf('All Trials (p=%0.1e)',p))
aH(2) = subplot(1,3,2,'parent',fH);
boxplot(aH(2),mFP_box12,{'box1 (-.)','box2 (--)'})
hold(aH(2),'on');
for k=1:size(mFP_box12,1)
    x1 = (rand-0.5)/3+1;
    x2 = (rand-0.5)/3+2;
    plot(aH(2),x1,mFP_box12(k,1),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);
    plot(aH(2),x2,mFP_box12(k,2),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);
    plot(aH(2),[x1,x2],mFP_box12(k,:),'LineStyle',':','Color',mrk_color{ptnum(k)})
end
[h,p,ci,stats] = ttest2(mFP_box12(:,1),mFP_box12(:,2));
ylabel(aH(2),'Normalized Aperiodic Exponent (Fooof)')
title(aH(2),sprintf('Mean Across Patient/Region (p=%0.1e)\nPatient (r=1, g=2, b=3, c=4, m=5)\nRegion (o=anthip, +=lattemp, *=entperi, x=posthip)',p))

aH(3) = subplot(1,3,3,'parent',fH);
dxx = d_ex./mean(d_ex,1);
dxx_me = mean(dxx,2,"omitnan");
dxx_se = 1.96*std(dxx,0,2,'omitnan')./sqrt(size(dxx,2)); %95ci
plot(aH(3),tsec_ex,dxx_me,'k')
hold on
plot(aH(3),tsec_ex,dxx_me+dxx_se,':k')
plot(aH(3),tsec_ex,dxx_me-dxx_se,':k')


%plotting full spectrum fit
power_box1 = squeeze(mean(FR_powerspec(tidx_box1_ex,:,:,plottype),1));
power_box2 = squeeze(mean(FR_powerspec(tidx_box2_ex,:,:,plottype),1));
fooofed_box1 = squeeze(mean(FR_fooofedspec(tidx_box1_ex,:,:,plottype),1));
fooofed_box2 = squeeze(mean(FR_fooofedspec(tidx_box2_ex,:,:,plottype),1));
apfit_box1 = squeeze(mean(FR_apfit(tidx_box1_ex,:,:,plottype),1));
apfit_box2 = squeeze(mean(FR_apfit(tidx_box2_ex,:,:,plottype),1));
r2_box1 = squeeze(mean(FR(tidx_box1_ex,6,:,plottype),1));
r2_box2 = squeeze(mean(FR(tidx_box2_ex,6,:,plottype),1));
err_box1 = squeeze(mean(FR(tidx_box1_ex,7,:,plottype),1));
err_box2 = squeeze(mean(FR(tidx_box2_ex,7,:,plottype),1));

for k=1:size(uPtCh,1)
    idx = (uPtChIdx==k);

    p_box1 = mean(power_box1(:,idx),2);
    p_box2 = mean(power_box2(:,idx),2);
    f_box1 = mean(fooofed_box1(:,idx),2);
    f_box2 = mean(fooofed_box2(:,idx),2);
    a_box1 = mean(apfit_box1(:,idx),2);
    a_box2 = mean(apfit_box2(:,idx),2);
    r_box1 = mean(r2_box1(idx));
    r_box2 = mean(r2_box2(idx));
    e_box1 = mean(err_box1(idx));
    e_box2 = mean(err_box2(idx));

    figure('Position',[50,50,1300,500]);
    subplot(1,2,1)
    plot(freqs,p_box1);
    hold on;
    plot(freqs,f_box1)
    plot(freqs,a_box1)
    xlim([freqs(1),freqs(end)])
    % set(gca,'XScale','linear')
    set(gca,'XScale','log')
    title(sprintf('Box1 (Pt=%d, R2=%0.2f, Err=%0.2f)',k,r_box1,e_box1))
    legend({'power','foooffit','apfit'})

    subplot(1,2,2)
    plot(freqs,p_box2);
    hold on;
    plot(freqs,f_box2)
    plot(freqs,a_box2)
    xlim([freqs(1),freqs(end)])
    % set(gca,'XScale','linear')
    set(gca,'XScale','log')
    title(sprintf('Box2 (Pt=%d, R2=%0.2f, Err=%0.2f)',k,r_box2,e_box2))
    legend({'power','foooffit','apfit'})

end


%Comparing fixed/knee R2
% R2_fixed = reshape(FR(:,6,:,1),[],1);
% R2_knee = reshape(FR(:,6,:,2),[],1);
% 
% [p,tbl,stats] = kruskalwallis([R2_fixed,R2_knee],{'fixed','knee'},'off'); 
% figure;
% [c,m] = multcompare(stats,'alpha',0.05,'display','on');

%Comparing fixed/knee variance of residuals
fidx = freqs<60;
Rv_fixed = reshape(var(FR_powerspec(:,fidx,:,1)-FR_fooofedspec(:,fidx,:,1),0,2),[],1);
Rv_knee = reshape(var(FR_powerspec(:,fidx,:,2)-FR_fooofedspec(:,fidx,:,2),0,2),[],1);

[p,tbl,stats] = kruskalwallis([Rv_fixed,Rv_knee],{'fixed','knee'},'off'); %use raw responses here since these are stats by subject
figure('position',[50,50,1600,650]);
subplot(1,2,1)
[c,m] = multcompare(stats,'alpha',0.05,'display','on');
title(sprintf('Variance of residuals (%0.0f to %0.0fHz)',min(freqs(fidx)),max(freqs(fidx))))

%Comparing fixed/knee of calculated R2 (these are lower than the values provided by fooof but show similar patterns and differences between fixed/knee)
SS_tot = reshape(sum((FR_powerspec(:,fidx,:,1)-mean(FR_powerspec(:,fidx,:,1),2)).^2,2),[],1);
SS_res = reshape(sum((FR_powerspec(:,fidx,:,1)-FR_fooofedspec(:,fidx,:,1)).^2,2),[],1);
R2_fixed_calc = 1-SS_res./SS_tot;

SS_tot = reshape(sum((FR_powerspec(:,fidx,:,2)-mean(FR_powerspec(:,fidx,:,2),2)).^2,2),[],1);
SS_res = reshape(sum((FR_powerspec(:,fidx,:,2)-FR_fooofedspec(:,fidx,:,2)).^2,2),[],1);
R2_knee_calc = 1-SS_res./SS_tot;

% figure;
% plot(R2_fixed);
% hold on;
% plot(R2_fixed_calc)

[p,tbl,stats] = kruskalwallis([R2_fixed_calc,R2_knee_calc],{'fixed','knee'},'off'); %use raw responses here since these are stats by subject
subplot(1,2,2)
[c,m] = multcompare(stats,'alpha',0.05,'display','on');
title(sprintf('R2 (%0.0f to %0.0fHz)',min(freqs(fidx)),max(freqs(fidx))))

