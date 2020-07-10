% The optimal algorithm reached after analysis of the 17 records was:
% D = .025 x [(.15T(i - 4) + .15T(i - 3) + .15T(i - 2) + .08T(i - 1) + .21T(i) + .12T(i + 1) + .13T(i + 2)]
% where T(i) represents the maximal epoch value in minute i, etc. If D >= 1.0, the minute was scored wake;
% otherwise it was scored sleep.
function T = detect_sleepWake(T,nFilt)
nR = 5:numel(T.odba_max)-2;
maxFilt = medfilt1(T.odba_max,nFilt);
im = zeros(numel(nR),7);
im(:,1) = 0.15 * maxFilt(nR-4);
im(:,2) = 0.15 * maxFilt(nR-3);
im(:,3) = 0.15 * maxFilt(nR-2);
im(:,4) = 0.08 * maxFilt(nR-1);
im(:,5) = 0.21 * maxFilt(nR);
im(:,6) = 0.21 * maxFilt(nR+1);
im(:,7) = 0.21 * maxFilt(nR+2);

W = prctile(T.odba,95) * sum(im,2); % multiplier determines threshold
Db = zeros(numel(T.odba_max),1);
Db(W >= 1) = 1;
Db = circshift(Db,4);
Db(1:4) = Db(5);
Db(end:end-1) = Db(end-2);
T.awake = Db;

% % % % close all
% % % % ff(1200,800);
% % % % showAmt = 3600*4;
% % % % xlims = [1,numel(T.odba);10000,10000+showAmt;50000,50000+showAmt];
% % % % titles = {'Several Days','Night Snippet','Day Snippet'};
% % % % colors = lines(3);
% % % % for ii = 1:4
% % % %     if ii == 1
% % % %         subplot(2,2,[1,2]);
% % % %         xlimVal = xlims(ii,:);
% % % %         useTitle = titles{1};
% % % %     elseif ii == 2
% % % %         continue;
% % % %     else
% % % %         subplot(2,2,ii);
% % % %         xlimVal = xlims(ii-1,:);
% % % %         useTitle = titles{ii-1};
% % % %     end
% % % %     plot(T.odba,'-','color',[0 0 0 0.1]);
% % % %     hold on;
% % % %     plot(T.odba_max,'k-');
% % % %     ylabel('ODBA');
% % % %     ylim([0 16]);
% % % %     yyaxis right;
% % % %     plot(W,'-','color',colors(3,:));
% % % %     hold on;
% % % %     plot(Db,'-','color','r');
% % % %     ylabel('Webster classifier (W)');
% % % %     xlim(xlimVal);
% % % %     ylim([0 10]);
% % % %     xlabel('time (s)');
% % % %     legend('odba','mov max','raw W','bin W');
% % % %     title(useTitle);
% % % % end