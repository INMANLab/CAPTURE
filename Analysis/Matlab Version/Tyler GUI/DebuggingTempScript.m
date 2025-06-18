%%
    if isempty(p.permtype)
        error('permtype must be specified!');
    end

    warning('off','MATLAB:contour:ConstantData');
    
    transtype = p.transtype; %'Outdoor Beg', 'Outdoor End', 'Doorway'
    regiontype = p.regiontype; %'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'
    walktype = p.walktype; %'First Walks','Last Walks','Stop Walks','Go Walks','All Walks', '1,2,5'
    veltype = p.veltype; %'', 'High Change', 'Low Change'
    patienttype = p.patienttype; %'All Patients', '2,4'
    desctype = p.desctype; %close, open, etc. (optional, will skip if empty)
    permtype = p.permtype; %standard, zscore
    correctiontype = p.correctiontype; %cluster, pixel, fdr, fdr+cluster (if empty, cluster correction is skipped)
    if strcmp(correctiontype,'fdr+cluster') && strcmp(permtype,'standard')
        error('fdr+cluster cannot be performed with the standard permtype!');
    end
    transrng = p.transrng; %default [-10,10]
    normrng = p.normrng; %default [-10,10] must be within transrng
    if ~p.fullwalknorm
        if normrng(1)<transrng(1) || normrng(2)>transrng(2) || normrng(2)<=normrng(1)
            error('normrng is incorrect!');
        end
    end
    pval = p.pval; %default p=0.05
    pvalclust = p.pvalclust; %default 0.01
    plottrials = p.plottrials;
    trialsfreqrng = p.trialsfreqrng;
    plotpowercomp = p.plotpowercomp;
    plotfooofcomp = p.plotfooofcomp;
    boxcomprng = p.boxcomprng; %2x4 matrix where rows are each box and cols are low/high ranges for time/freq
    if isempty(p.clim)
        switch permtype
            case 'zscore'
                climit = [-10,10];
            case 'standard' %dB
                climit = [-1,1];
        end
    else
        climit = p.clim;
    end

    %Getting trials and specgram data
    obj.filterMultTransData(varargin{:}); %table of all trials and corresponding info (this creates MT and must be run before getFilteredMultTransData)
    obj.getFilteredMultTransData(transrng); %raw power for all trials (time x freq x trial)

    %Init some params
    tsamp = obj.MultTrans.MTd.tsamp_wv; %+-10sec at 25Hz (downsampled by 10 from 250)
    tsec = obj.MultTrans.MTd.tsec_wv;
    baseidx = tsec>=normrng(1) & tsec<=normrng(2); %baseline (normalization) indices
    nperm = 1000; %permutations

    pwr = obj.MultTrans.MTd.d_wv;
    tsec = obj.MultTrans.MTd.tsec_wv;
    freq = obj.MultTrans.MTd.freq_wv;
    np = obj.MultTrans.MTd.d_np;
    tsec_np = obj.MultTrans.MTd.tsec_np;

    predlist = [{'xs','gz','kd','am','ed','pe1','pe2','pe3','gy'};{'Vel','Fix','KDE','Amb','Eda','PeT','PeA','PeB','HeadTurn'}];

    nfreq = length(freq);
    ntime = length(tsamp);
    ntrials = size(obj.MultTrans.MT,1);

    %Real
    if p.fullwalknorm %Normalization is done across the entire walk for each chan in getMultTransData, so skip normalization across window
        bpwr = ones(1,nfreq); %Since the entire walk is baseline normalized, set to ones for full window norm to preserve relative power
    else
        bpwr = squeeze(mean(mean(pwr(baseidx,:,:),1,"omitnan"),3,"omitnan")); %mean across time and then trials (1 x freq)
        pwr = pwr./mean(pwr,1,"omitnan"); %full epoch norm (time x freq x trial) 
    end
    mpwr = squeeze(mean(pwr,3,"omitnan")); %trial avg (time x freq)
    npwr = 10*log10(mpwr./bpwr); %normalized power in dB
    npwr_trials = 10*log10(pwr./bpwr); %normalized power for individual trials
%
            %Permutations
            disp('Running permutations...')
            % PM = calcRWAPerm(pwr(baseidx,:,:),bpwr,0); %trialtype=0 for single specgram (vector of -1 or 1 for all trials when computing specgram difference)
            PM = calcRWAPerm_mex(pwr(baseidx,:,:),bpwr,0); %1.3 times faster

            %Finding percentiles of permuted distributions
            pm_freq = reshape(permute(PM,[1,3,2]),[],nfreq); %(ntime*nperm) x nfreq
            pc = prctile(pm_freq,[pval*100/2,100-pval*100/2],1); %2 x nfreq
            mPM = mean(pm_freq,1); %mean across permutations and time (1 x nfreq)
            sPM = std(pm_freq,0,1); %std across permutations and time (1 x nfreq)
%             %Keep time (similar to the zscore method)
%             pm = reshape(PM,[],nperm); %finding percentiles for all time x freq distributions (instead of collapsing across time like above)
%             pc_low = reshape(prctile(pm,pval*100/2,2),ntime,nfreq); %two tail (across 2nd dim -> 1000 perms)
%             pc_high = reshape(prctile(pm,100-pval*100/2,2),ntime,nfreq); %two tail
%             pc = permute(cat(3,pc_low,pc_high),[3,1,2]); %2 x ntime x nfreq
%             mPM = mean(PM,3); %mean across permutations
%             sPM = std(PM,0,3); %std across permutations

            %Finding clusters/max pixels in permuted data
            max_clust_info = nan(nperm,1);
            max_pixel_pvals = zeros(nperm,2);
            for m=1:nperm
                pm = PM(:,:,m);
                switch permtype
                    case 'standard'
                        max_pixel_pvals(m,:) = [min(pm(:)),max(pm(:))]; %pixel correction distributions (pooled across time and freq)
                        pm(pm>squeeze(pc(1,:,:)) & pm<squeeze(pc(2,:,:))) = 0;
                    case 'zscore'
                        pm = (pm-mPM)./sPM;
                        max_pixel_pvals(m,:) = [min(pm(:)),max(pm(:))]; %pixel correction distributions (pooled across time and freq)
                        switch correctiontype
                            case 'fdr+cluster'
                                PV = 1-normcdf(abs(pm));
                                PV = reshape(mafdr(PV(:),'BHFDR',true),size(PV));
                                pm(PV>=pval) = 0;
                            otherwise
                                pm(abs(pm)<norminv(1-pval)) = 0;
                        end
                end
                clustinfo = bwconncomp(pm);
                stats = regionprops(clustinfo,pm,'pixelvalues');
                max_clust_info(m) = max([0,abs(cellfun(@sum,{stats.PixelValues}))]); %max cluster size using pixel values for each permutation
                % max_clust_info(m) = max([0,cellfun(@numel,clustinfo.PixelIdxList)]); %max cluster size for each permutation
            end

            %Finding percentiles of pixel-level distributions
            pc_pixel = prctile(max_pixel_pvals(:),[pval*100/2,100-pval*100/2]); %two tail

            %Thresholding real
            switch permtype
                case 'standard'
                    npwr_thresh = npwr;
                    switch correctiontype
                        case 'fdr' %results are slightly different than PermSpecGram because this operates on log10 data and does rotating shuffle
                            %Finding pvals based on rank in distribution for FDR correction
                            % PV = calcRWAPVal(npwr,pm_freq);
                            PV = calcRWAPVal_mex(npwr,pm_freq); %~3x faster
                            PV = reshape(mafdr(PV(:),'BHFDR',true),size(npwr));
                            npwr_thresh(PV>=pval) = 0;
                        case 'pixel'
                            npwr_thresh(npwr_thresh>pc_pixel(1) & npwr_thresh<pc_pixel(2)) = 0;
                        otherwise %cluster or empty
                            npwr_thresh(npwr_thresh>squeeze(pc(1,:,:)) & npwr_thresh<squeeze(pc(2,:,:))) = 0;
                    end
                case 'zscore'
                    npwr = (npwr-mPM)./sPM;
                    npwr_thresh = npwr;
                    npwr_trials = (npwr_trials-mPM)./sPM;
                    switch correctiontype
                        case {'fdr','fdr+cluster'}
                            PV = 1-normcdf(abs(npwr));
                            PV = reshape(mafdr(PV(:),'BHFDR',true),size(PV));
                            npwr_thresh(PV>=pval) = 0;
                        case 'pixel'
                            npwr_thresh(npwr_thresh>pc_pixel(1) & npwr_thresh<pc_pixel(2)) = 0;
                        otherwise %cluster or empty
                            npwr_thresh(abs(npwr)<norminv(1-pval))=0; %norminv gives stat value at pval for normal distribution (i.e. p=0.05 is -1.65)
                    end
            end

            %Removing clusters
            if contains(correctiontype,'cluster') && ~isempty(correctiontype)
                if size(PM,1)~=ntime
                    error('Cluster correction cannot be performed if the baseline/normalization window is smaller than the transition window!');
                else
                    clustinfo = bwconncomp(npwr_thresh);
                    stats = regionprops(clustinfo,npwr_thresh,'pixelvalues');
                    clust_info = abs(cellfun(@sum,{stats.PixelValues}));
                    % clust_info = cellfun(@numel,clustinfo.PixelIdxList);
                    clust_threshold = prctile(max_clust_info,100-pvalclust*100);
                    whichclusters2remove = find(clust_info<clust_threshold);
                    for k=1:length(whichclusters2remove)
                        npwr_thresh(clustinfo.PixelIdxList{whichclusters2remove(k)})=0;
                    end
                end
            end

            %Plotting
            fH = figure('Position',[50,50,1200,800]);
            aH = axes('parent',fH);
            colormap("jet");
            contourf(tsec,freq,npwr',100,'linecolor','none','parent',aH);
            hold(aH,"on");
            contour(tsec,freq,(npwr_thresh~=0)',1,'parent',aH,'linecolor','w','linewidth',2);
            set(aH,'yscale','log','YTick',2.^(1:6),'yticklabel',2.^(1:6),'yminortick','off','clim',climit); %now in units of standard deviation
            plot(aH,[0,0],[2,120],'--k','LineWidth',2);
            pH = nan(size(predlist,2),3);
            yyaxis(aH,'right');
            for k=1:size(predlist,2) %xs,gz,kd,am,ed,pe1,pe2,pe3,gy
                dx = obj.MultTrans.MTd.(['d_',predlist{1,k}]);
                tx = obj.MultTrans.MTd.(['tsec_',predlist{1,k}]);
                if contains(predlist{1,k},'pe') || contains(predlist{1,k},'gy')
                    dxx = dx;
                else
                    dxx = dx(:,~isoutlier(mean(dx)));
                end
                dxx_me = mean(dxx,2,"omitnan");
                dxx_se = 1.96*std(dxx,0,2,'omitnan')./sqrt(size(dxx,2)); %95ci
                pH(k,1) = plot(aH,tx,dxx_me,'-k','LineWidth',2);
                pH(k,2) = plot(aH,tx,dxx_me-dxx_se,':k','LineWidth',0.5);
                pH(k,3) = plot(aH,tx,dxx_me+dxx_se,':k','LineWidth',0.5);
            end  
            predidx = strcmp(predlist(2,:),p.predtype);
            set(pH(~predidx,:),'visible','off');
            ylabel(aH,p.predtype);
            switch p.predtype %vel:0-1.5, gaze:0-1, kde:1-7, amb: 0-500, eda:-0.2-0.5, pe:0-1
                case 'Vel'
                    ylim(aH,[0,1.5])
                case 'Fix'
                    ylim(aH,[0,1])
                case 'KDE'
                    ylim(aH,[1,7])
                case 'Amb'
                    ylim(aH,[0,500])
                case 'Eda'
                    ylim(aH,[-0.2,0.5])
                case {'PeT','PeA','PeB'}
                    ylim(aH,[0,0.3])
                case 'HeadTurn' 
                    ylim(aH,[-40,40])
                otherwise
                    aH.YAxis(2).Visible = 'off';
            end
            yyaxis(aH,'left');
            if plottrials
                plot(aH,[tsec(1),tsec(end)],repmat(trialsfreqrng,2,1),':k')
            end
            if plotpowercomp || plotfooofcomp
                rectangle(aH,'Position',[boxcomprng(1,1),boxcomprng(1,3),diff(boxcomprng(1,1:2)),diff(boxcomprng(1,3:4))],'LineStyle','-.');
                rectangle(aH,'Position',[boxcomprng(2,1),boxcomprng(2,3),diff(boxcomprng(2,1:2)),diff(boxcomprng(2,3:4))],'LineStyle','--');
            end
            cb = colorbar(aH);
            cblims = [cb.Limits(1),0,cb.Limits(2)];
            cb.Ticks = cblims;
            cb.TickLabels = cblims;
            xlim(aH,[tsec(1),tsec(end)])
            ylim(aH,[2,120])
            if contains(permtype,'zscore')
                ylabel(cb,'Zscore')
            else
                ylabel(cb,'dB')
            end
            xlabel(aH,'sec');
            ylabel(aH,'Hz');
            if isempty(str2num(patienttype))
                pttypestr = patienttype;
            else
                pttypestr = 'Sel Patients';
            end
            ptstr = regexprep(num2str(unique([obj.MultTrans.MT.patient])'),'\s+','');
            if isempty(str2num(walktype))
                wktypestr = walktype;
            else
                wktypestr = 'Sel Walks';
            end
            wkstr = regexprep(num2str(unique([obj.MultTrans.MT.walk])'),'\s+','');
            rgcell = {'a','l','e','p'}; %AntHipp, LatTemp, Ent+Peri, PostHipp+Para
            rgstr = cell2mat(rgcell(unique(obj.MultTrans.MT.regionnum)));
            if p.fullwalknorm
                normstr = 'full';
            else
                normstr = regexprep(num2str(normrng),'\s+','t');
            end
            transstr = regexprep(num2str(transrng),'\s+','t');
            [PercOverlap,AvgOverlapSec] = obj.calcOverlap(obj.MultTrans.MT,transrng);
            ttlstr = sprintf('%s, %s, %s-%s, %s-%s, %s-%s, %s \n(n=%s, p<%1.2f, %s, %s, b=%s, w=%s, o=%0.0fp-%0.0fs)',...
                transtype,desctype,regiontype,rgstr,pttypestr,ptstr,...
                wktypestr,wkstr,veltype,num2str(ntrials),...
                pval,permtype,correctiontype,normstr,transstr,PercOverlap,AvgOverlapSec);
            title(aH,ttlstr); 

            fH2 = [];
            if plottrials
                fH2 = figure('Position',[50,50,1200,800]);
                aH2 = axes('parent',fH2);
                fH2.Colormap = colormap('jet');

                [mt,sidx] = sortrows(obj.MultTrans.MT,{'patient','chan','walk'});
                fidx = freq>=trialsfreqrng(1) & freq<=trialsfreqrng(2);
                fdat = squeeze(mean(npwr_trials(:,fidx,sidx),2))';

                cl = round(prctile(abs(fdat(:)),99));
                CLimit2 = [-cl,cl];
                II = fix((fdat-CLimit2(1))./(CLimit2(2)-CLimit2(1))*255)+1;
                II(II<1) = 1; II(II>256) = 256;
                iH = imagesc(tsec,1:size(fdat,1),II,'parent',aH2);
                iH.CDataMapping = "direct";

                xlim(aH2,[tsec(1),tsec(end)])
                cb2 = colorbar(aH2);
                cb2.Limits = [cb2.Limits(1),cb2.Limits(2)-1];
                cblims2 = [cb2.Limits(1),(diff(cb2.Limits)-1)/2,cb2.Limits(2)-1];
                cb2.Ticks = cblims2;
                cb2.TickLabels = [CLimit2(1),0,CLimit2(2)];
                if contains(permtype,'zscore')
                    ylabel(cb2,'Zscore')
                else
                    ylabel(cb2,'dB')
                end
                title(aH2,sprintf('Trials %0.0f-%0.0f Hz (%s)',trialsfreqrng(1),trialsfreqrng(2),ttlstr))
                [uPt,~,uPtIdx] = unique(mt.patient);
                mYtick = nan(length(uPt),1);
                hold(aH2,'on');
                for k=1:length(uPt)
                    idx = find(uPtIdx==k);
                    mYtick(k) = mean(idx);
                    plot([tsec(1)+0.3,tsec(1)+0.3],[idx(1),idx(end)],'color',[0.5,0.5,0.5],'LineWidth',10)
                end
                set(aH2,'ytick',mYtick,'yticklabel',uPt)
                ylabel(aH2,'patient')
                xlabel(aH2,'sec')
            end

            [uPtCh,~,uPtChIdx] = unique(obj.MultTrans.MT(:,{'patient','regionnum'}));
            mrk_style = {'o','+','*','x'}; %by region -> 1=AntHipp, 2=LatTemp, 3=Ent+Peri, 4=PostHipp+Para
            mrk_size = [6,10,10,10];
            mrk_color = {'r','g','b','c','m'}; %by patient -> 1:5
            ptnum = uPtCh.patient;
            rgnum = uPtCh.regionnum;

            fH3 = [];
            if plotpowercomp
                tidx_box1 = tsec>boxcomprng(1,1) & tsec<boxcomprng(1,2);
                tidx_box2 = tsec>boxcomprng(2,1) & tsec<boxcomprng(2,2);
                fidx_box1 = freq>boxcomprng(1,3) & freq<boxcomprng(1,4);
                fidx_box2 = freq>boxcomprng(2,3) & freq<boxcomprng(2,4);

                % pwr_box1 = squeeze(mean(mean(10*log10(pwr(tidx_box1,fidx_box1,:)),1),2));
                % pwr_box2 = squeeze(mean(mean(10*log10(pwr(tidx_box2,fidx_box2,:)),1),2));

                pwr_box1 = squeeze(mean(mean(pwr(tidx_box1,fidx_box1,:),1),2)); %average before converting to dB
                pwr_box2 = squeeze(mean(mean(pwr(tidx_box2,fidx_box2,:),1),2));

                % tidx_xs_box1 = tsec_xs>boxcomprng(1,1) & tsec_xs<boxcomprng(1,2);
                % tidx_xs_box2 = tsec_xs>boxcomprng(2,1) & tsec_xs<boxcomprng(2,2);

                % dxx = dx;
                % dxx(:,isoutlier(mean(dx))) = nan;
                % dx_box1 = mean(dxx(tidx_xs_box1,:))';
                % dx_box2 = mean(dxx(tidx_xs_box2,:))';

                mpwr_box1 = nan(size(uPtCh,1),1);
                mpwr_box2 = nan(size(uPtCh,1),1);
                % ccdx_box1 = nan(size(uPtCh,1),1);
                % ccdx_box2 = nan(size(uPtCh,1),1);
                for k=1:size(uPtCh,1)
                    idx = (uPtChIdx==k);
                    pwr_b1 = pwr_box1(idx);
                    pwr_b2 = pwr_box2(idx);
                    % dx_b1 = dx_box1(idx);
                    % dx_b2 = dx_box2(idx);
                    mpwr_box1(k) = mean(pwr_b1);
                    mpwr_box2(k) = mean(pwr_b2);
                    % nan_idx = isnan(dx_b1)|isnan(dx_b2);
                    % cc_b1 = corrcoef(dx_b1(~nan_idx),10*log10(pwr_b1(~nan_idx))); %convert pwr to dB before correlation
                    % cc_b2 = corrcoef(dx_b2(~nan_idx),10*log10(pwr_b2(~nan_idx)));
                    % ccdx_box1(k) = cc_b1(1,2);
                    % ccdx_box2(k) = cc_b2(1,2);
                end

                mpwr_box12 = median([pwr_box1,pwr_box2]);
                pwr_chg_box12 = (mpwr_box12(2)-mpwr_box12(1))/mpwr_box12(1); %percentage change in normalized power

                mmpwr_box12 = median([mpwr_box1,mpwr_box2]);
                mpwr_chg_box12 = (mmpwr_box12(2)-mmpwr_box12(1))/mmpwr_box12(1); %percentage change in normalized power

                fH3 = figure('Position',[50,50,1200,550]); 
                aH3 = [];
                aH3(1) = subplot(1,2,1,'parent',fH3);
                boxplot(aH3(1),10*log10([pwr_box1,pwr_box2]),{'box1 (-.)','box2 (--)'})
                [h,p,ci,stats] = ttest2(10*log10(pwr_box1),10*log10(pwr_box2));
                ylabel(aH3(1),'Power (dB)')
                title(aH3(1),sprintf('All Trials (p=%0.1e, %0.0f%%)',p,pwr_chg_box12*100))
                aH3(2) = subplot(1,2,2,'parent',fH3);
                boxplot(aH3(2),10*log10([mpwr_box1,mpwr_box2]),{'box1 (-.)','box2 (--)'})
                hold(aH3(2),'on');
                for k=1:length(mpwr_box1)
                    x1 = (rand-0.5)/3+1;
                    x2 = (rand-0.5)/3+2;
                    plot(aH3(2),x1,10*log10(mpwr_box1(k)),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);
                    plot(aH3(2),x2,10*log10(mpwr_box2(k)),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);                    
                    plot(aH3(2),[x1,x2],10*log10([mpwr_box1(k),mpwr_box2(k)]),'LineStyle',':','Color',mrk_color{ptnum(k)})
                end
                [h,p,ci,stats] = ttest2(10*log10(mpwr_box1),10*log10(mpwr_box2));
                ylabel(aH3(2),'Mean Power (dB)')
                title(aH3(2),sprintf('Mean Across Patient/Region (p=%0.1e, %0.0f%%)\nPatient (r=1, g=2, b=3, c=4, m=5)\nRegion (o=anthip, +=lattemp, *=entperi, x=posthip)',p,mpwr_chg_box12*100))
            end %power comparison
%%
            fH4 = [];

            

                


%%

