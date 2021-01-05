% Fpass = 0.7;
% Fstop = 0.5;
% Apass = 1;
% Astop = 65;
% Fs = 200;
% 
% d = designfilt('highpassfir', ...
%   'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
%   'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
%   'DesignMethod','equiripple','SampleRate',Fs);
% 
% fvtool(d)

%% �ο�matlab�����źŴ��� 85�� ������
% ���ȥ�������������С���˷���ϣ���ϳ��Ķ���ʽΪ������
clc; clear; close all;
load ECG.mat
plot(BCG)
BCG1 = BCG1 - mean(BCG1);
BCG1 = BCG1(500:2999);

Fs = 200;
N = length(BCG1);
t = (0:N-1)'/Fs;

a = polyfit(t, BCG1, 3);
xtrend = polyval(a, t);

y1 = sgolayfilt(BCG1, 5, 501);

subplot 311;
plot(BCG1 - xtrend);

subplot 312;
plot(BCG1);

subplot 313;
plot(detrend(BCG1));

%% median filter
% % [1]��ΰ��,�봺.һ��ʵ�õ�ȥ����Ư���˲��㷨[J].���ݴ�ѧѧ��(���ư�),2006(01):62-64.
% clc; clear; close all;
% load ECG.mat
% 
% BCG = BCG';
% BCG_max = max(BCG);
% BCG_min = min(BCG);
% BCG = ((BCG - BCG_min) / (BCG_max - BCG_min)) * 2 - 1;
% 
% 
% W = 51;
% L = length(BCG);
% 
% BCG_expand(1:(W-1)/2) = BCG(1);
% BCG_expand((W-1)/2+1:L+(W-1)/2) = BCG(1:L);
% BCG_expand(L+(W-1)/2+1:L+(W-1)) = BCG(L);
% 
% BCG_baseline(1:L) = 0;
% for i = 1:1:L
%     BCG_baseline(i) = median(BCG_expand(i:i+W-1));
% end
% 
% BCG_filtered = BCG - BCG_baseline;
% 
% plot(BCG_baseline); hold on;
% plot(BCG); hold on;
% plot(BCG_filtered);

%% cite other people's work
% clc; clear; close all;
% load ECG.mat
% 
% BCG_max = max(BCG);
% BCG_min = min(BCG);
% BCG = ((BCG - BCG_min) / (BCG_max - BCG_min));
% BCG = BCG;
% Fs = 200;
% 
% % BLremover_instance = BLremover;
% % BCG_filtered = BLremover_instance.IIRRemoveBL(BCG, Fs);
% BCG_filtered = FARemoveBL(BCG, 200);
% plot(BCG_filtered); hold on;
% plot(BCG)


%% wavelet filter to remove baseline��
% % if baseline wander in large sacle, this method does not fit
% % [1]��С��,����,�����,ղ����.����С���任������Ӧ�˲�������ECG�л���Ư��[J].�й���ѧ������ѧѧ��,2000(04):75-79.
% % [1]ʷ����,�ƽ���,��Ӣ��,�ƽ���.����С��������Ӧ�˲���ECG����Ư��У��[J].���������,2013,39(11):226-229+244.
% clc; clear; close all;
% load ECG.mat
% BCG_max = max(BCG);
% BCG_min = min(BCG);
% BCG = ((BCG - BCG_min) / (BCG_max - BCG_min));
%  
% N = length(BCG);
% leval = 7;
% 
% [C, L] = wavedec(BCG, leval, 'db1');
% % Ca = appcoef(c,l,'db1', 1);
% % Cd = detcoef(c,l,[1:1:9]);
% % A1 = upcoef('a',Ca,'db1',1,N);
% A10 = wrcoef('a',C,L,'db1',leval);
% 
% for i = 1:1:leval
%     D(:, i) = wrcoef('d',C,L,'db1',i);
% end
% 
% W = 0.01*ones(size(D,2),1);
% lr = 0.001;
% result = BCG;
% loss = 0.5;
% i = 0;
% 
% while(i<1000 && loss<5)
%     i = i+1;
%     error = (result - D*W);
%     W = W + 2 * lr * D' * error;
%     
%     loss = sum((result - D*W).^2);
% end
% 
% plot(D*W)
% 
% % result = zeros(N, 1);
% % for i=1:1:leval
% %     result = result + D(:, i);
% % end
% 
% % plot(result)






