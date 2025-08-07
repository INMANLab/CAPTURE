function FooofPlotGif(np, tsec_np,boxcomprng) 
% FoooF Peak falls in which EEG band.
fRange1 = boxcomprng(1,3):.2:boxcomprng(1,4);
fRange2 = boxcomprng(2,3):.2:boxcomprng(2,4);

tidx1 = tsec_np>boxcomprng(1,1) & tsec_np<boxcomprng(1,2);
tidx2 = tsec_np>boxcomprng(2,1) & tsec_np<boxcomprng(2,2);

settings.max_n_peaks = 3; %aperiodic component has a better fit if more peaks are removed
settings.aperiodic_mode = 'knee';
settings.aperiodic_mode = 'fixed';


signal1 = np(tidx1,:);
signal2 = np(tidx2,:);
psd1 = pwelch(signal1, 100, [], fRange1, 250);
psd2 = pwelch(signal1, 100, [], fRange2, 250);

pWave1 = zeros(length(fRange1),size(np,2));
pWave2 = zeros(length(fRange2),size(np,2));
for pIdx = 1:size(np,2)
    pWave1(:,pIdx) = mean(abs(calcWavTF(signal1(:,pIdx),fRange1,250)).^2,2);
    pWave1(:,pIdx) = pWave1(:,pIdx)/mean(pWave1(:,pIdx));
    pWave2(:,pIdx) = mean(abs(calcWavTF(signal2(:,pIdx),fRange2,250)).^2,2);
    pWave2(:,pIdx) = pWave2(:,pIdx)/mean(pWave2(:,pIdx));
end

pAv_Welch1 = mean(psd1,2);
pAv_Welch2 = mean(psd2,2);
pAv_Wave1  = mean(pWave1,2);
pAv_Wave2  = mean(pWave2,2);


figure
subplot 121
hold on
plot(fRange1,10*log10(pAv_Welch1))
plot(fRange1,10*log10(pAv_Wave1))
legend("Wavelet","Welch")
xlabel("Freq")
title("First Half")

subplot 122
hold on
plot(log10(fRange1),10*log10(pAv_Welch1))
plot(log10(fRange1),10*log10(pAv_Wave1))
legend("Wavelet","Welch")
xlabel("log(Freq)")
title("first Half")

figure
subplot 121
hold on
plot(fRange2,10*log10(pAv_Welch2))
plot(fRange2,10*log10(pAv_Wave2))
legend("Wavelet","Welch")
xlabel("Freq")
title("Second Half")

subplot 122
hold on
plot(log10(fRange2),10*log10(pAv_Welch2))
plot(log10(fRange2),10*log10(pAv_Wave2))
legend("Wavelet","Welch")
xlabel("log(Freq)")
title("Second Half")


fooof_results2 = fooof(fRange1, pAv_Wave1', [fRange1(1), fRange1(end)], settings, true);
fooof_results2 = fooof(fRange2, pAv_Wave1', [fRange2(1), fRange2(end)], settings, true);

figure
subplot 121
fooof_plot(fooof_results1)
subplot 122
fooof_plot(fooof_results2)


a=0;

segNum = 1;
figure
subplot 121
hold on
plot(fRange1,10*log10(psd1(:,segNum)'))
plot(fRange1,10*log10(pWave1(:,segNum)'))
legend("Wavelet","Welch")
xlabel("Freq")
title("First Half")

subplot 122
hold on
plot(fRange1,10*log10(psd2(:,segNum)'))
plot(fRange1,10*log10(pWave2(:,segNum)'))
legend("Wavelet","Welch")
xlabel("log(Freq)")
title("first Half")




% 
% 
% 
% fbins = 2:0.2:85;    
% winSize= 500;   
% tWin = false(length(tsec_np),1);
% tWin(1:winSize) = true;
% 
% settings.max_n_peaks = 2; %aperiodic component has a better fit if more peaks are removed
% settings.aperiodic_mode = 'fixed';
% 
% v = VideoWriter("testVid", 'MPEG-4');
% v.FrameRate = 5;  % Adjust for desired speed
% open(v);
% 
% fig = figure('Visible', 'on');
% 
% for itt = 1:200:(length(tsec_np)-(winSize))
%     cla(fig)
%     [psda, freqsa] = pwelch(np(tWin,:), [], [], fbins, 250);
%     freqsa = freqsa';
%     psda = psda';
% 
%     psda = mean(psda);
% 
%     % Run FOOOF, also returning the model
%     fooof_results = fooof(freqsa, psda, [fbins(1), fbins(end)], settings, true);
%     fooof_plot(fooof_results)
% 
%     % pause(.5)
%      % Write frame to video
%     frame = getframe(fig);
%     writeVideo(v, frame);
% 
%     tWin = circshift(tWin,200);
% end
% 
% close(v);
% % disp("Video saved in "+ fileName);
% % implay('testVid.mp4')
% 
% tidx_box1_np = tsec_np>boxcomprng(1,1) & tsec_np<boxcomprng(1,2);
% tidx_box2_np = tsec_np>boxcomprng(2,1) & tsec_np<boxcomprng(2,2);
% 
% % fbins = 2:0.2:85;
% % pwr_box1 = mean(abs(calcWavTF(np(tidx_box1_np,1),fbins,250)).^2,2);
% % figure;plot(fbins,db(pwr_box1))
% % [psd1, freqs1] = pwelch(np(tidx_box1_np,1), [], [], fbins, 250);
% % figure;plot(freqs1,db(psd1))
% % pwr_box2 = mean(abs(calcWavTF(np(tidx_box2_np,:),fbins,250)).^2,2);
% 
% [psd1, freqs1] = pwelch(np(tidx_box1_np,:), [], [], fbins, 250);
% [psd2, freqs2] = pwelch(np(tidx_box2_np,:), [], [], fbins, 250);
% 
% freqs1 = freqs1';
% psd1 = psd1';
% freqs2 = freqs2';
% psd2 = psd2';
% psd1 = mean(psd1);
% psd2 = mean(psd2);
% 
% fooof_results1 = fooof(freqs1, psd1, [fbins(1), fbins(end)], settings, true);
% fooof_results2 = fooof(freqs2, psd2, [fbins(1), fbins(end)], settings, true);
% 
% figure
% subplot 121
% fooof_plot(fooof_results1)
% subplot 122
% fooof_plot(fooof_results2)
% 
% % 
% % FP_box1 = nan(size(np,2),length(FPLabels));
% % FP_box2 = nan(size(np,2),length(FPLabels));
% % wait_msg = parfor_wait(size(np,2));
% % for k=1:size(np,2)
% % for k=1:size(np,2)
% %     % wait_msg.Send;
% %     pwr_box1 = mean(abs(calcWavTF(np(tidx_box1_np,k),fbins,250)).^2,2);
% %     pwr_box2 = mean(abs(calcWavTF(np(tidx_box2_np,k),fbins,250)).^2,2);
% % 
% %     [psd, freqs] = pwelch(np(tidx_box1_np,k), 500, [], fbins, 250);
% %     [psda, freqsa] = pwelch(np(tidx_box1_np,:), 500, [], fbins, 250);
% % 
% %     % Transpose, to make inputs row vectors
% %     freqsa = freqsa';
% %     psda = psda';
% % 
% %     psda = mean(psda);
% % 
% %     % Run FOOOF, also returning the model
% %     fooof_results = fooof(freqsa, psda, [fbins(1), fbins(end)], settings, true);
% % 
% %     % Plot the resulting model
% %     fooof_plot(fooof_results)
% % 
% %     fr_box1 = fooof(fbins, pwr_box1, [fbins(1),fbins(end)], settings, 1); %fr.peak_params = [freq, height (aperiodic removed), width]
% %     fr_box2 = fooof(fbins, pwr_box2, [fbins(1),fbins(end)], settings, 1); %fr.peak_params = [freq, height (aperiodic removed), width]
% %     if isempty(fr_box1.peak_params)
% %         pkp = nan(1,3);
% %     else
% %         [~,midx] = max(fr_box1.peak_params(:,2)); %find the highest peak
% %         pkp = fr_box1.peak_params(midx,:);
% %     end
% %     FP_box1(k,:) = [pkp,fr_box1.aperiodic_params,fr_box1.r_squared,fr_box1.error];
% %     if isempty(fr_box2.peak_params)
% %         pkp = nan(1,3);
% %     else
% %         [~,midx] = max(fr_box2.peak_params(:,2)); %find the highest peak
% %         pkp = fr_box2.peak_params(midx,:);
% %     end
% %     FP_box2(k,:) = [pkp,fr_box2.aperiodic_params,fr_box2.r_squared,fr_box2.error];
% % end
% % % wait_msg.Destroy;
% % 
% % nFP_box12 = nan(size(FP_box1,1),2); %mean normalized across box1/box2 by patient/chan
% % mFP_box12 = nan(size(uPtCh,1),2);
% % for k=1:size(uPtCh,1)
% %     idx = (uPtChIdx==k);
% %     fp_box12 = [FP_box1(idx,5),FP_box2(idx,5)];
% %     nfp_box12 = fp_box12./mean(fp_box12,2);
% %     nFP_box12(idx,:) = nfp_box12;
% %     mFP_box12(k,:) = mean(nfp_box12,1);
% % end
% % 
% % % fH4 = figure('Position',[50,50,1200,550]);
% % % aH4 = [];
% % % aH4(1) = subplot(1,2,1,'parent',fH4);
% % % boxplot(aH4(1),nFP_box12,{'box1 (-.)','box2 (--)'})
% % % [h,p,ci,stats] = ttest2(nFP_box12(:,1),nFP_box12(:,2));
% % % ylabel(aH4(1),'Normalized Aperiodic Exponent (Fooof)')
% % % title(aH4(1),sprintf('All Trials (p=%0.1e)',p))
% % % aH4(2) = subplot(1,2,2,'parent',fH4);
% % % boxplot(aH4(2),mFP_box12,{'box1 (-.)','box2 (--)'})
% % % hold(aH4(2),'on');
% % % for k=1:size(mFP_box12,1)
% % %     x1 = (rand-0.5)/3+1;
% % %     x2 = (rand-0.5)/3+2;
% % %     plot(aH4(2),x1,mFP_box12(k,1),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);
% % %     plot(aH4(2),x2,mFP_box12(k,2),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);
% % %     plot(aH4(2),[x1,x2],mFP_box12(k,:),'LineStyle',':','Color',mrk_color{ptnum(k)})
% % % end
% % % [h,p,ci,stats] = ttest2(mFP_box12(:,1),mFP_box12(:,2));
% % % ylabel(aH4(2),'Normalized Aperiodic Exponent (Fooof)')
% % % title(aH4(2),sprintf('Mean Across Patient/Region (p=%0.1e)\nPatient (r=1, g=2, b=3, c=4, m=5)\nRegion (o=anthip, +=lattemp, *=entperi, x=posthip)',p))
% end
% 
% 
% %%
% % baseline = mean(pwr_box1);
% % pT = pwr_box1/baseline;
% % nWin = 100;
% % p1 = pwelch(np(tidx_box1_np,k),nWin,[],fbins,250);
% % p1 = p1/mean(p1);
% % 
% % figure
% % subplot 121
% % hold on
% % plot(fbins,10*log10(pT))
% % plot(fbins,10*log10(p1))
% % legend("Morlet","Welch")
% % 
% % 
% % subplot 122
% % hold on
% % plot(log(fbins),10*log10(pT))
% % plot(log(fbins),10*log10(p1))
% % legend("Morlet","Welch")
