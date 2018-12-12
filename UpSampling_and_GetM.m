function [data_Peak,SNR_tmp] = UpSampling_and_GetM(s_ac_test,Np)
% 用来对某定标器、某通道的数据切片进行二维升采样，然后取出峰值点的幅相信息
%
% 输入变量：
%   1）s_ac_test 是某定标器在某通道下的数据切片
%   2）Np 是该切片的大小：Np×Np
%
% 输出变量：
%   1）data_Peak 是该定标器的峰值点在该通道下的幅相信息
%   2）SNR_tmp   是该定标器在该通道下的 SNR 估计结果
%
% 该程序截止至：2016.10.31. 15:49

%%
% （1）进行二维升采样
Up_Sampling_Ratio = 16;                 % 二维升采样时，每一维的升采样倍数

S_ac_test_1 = fft(s_ac_test,[],1);      % 方位向fft
S_ac_test_2 = fft(S_ac_test_1,[],2);    % 距离向fft
clear s_ac_test;clear S_ac_test_1;

% 接下来进行二维补零
S_ac_test_buling_1 = zeros(Np,Up_Sampling_Ratio*Np);% 中间变量
S_ac_test_buling = zeros(Up_Sampling_Ratio*Np,Up_Sampling_Ratio*Np);
% ========================================================================
S_ac_test_buling_1(:,1:Np/2) = S_ac_test_2(:,1:Np/2);
S_ac_test_buling_1(:,Up_Sampling_Ratio*Np-(Np-Np/2)+1:Up_Sampling_Ratio*Np) =...
    S_ac_test_2(:,Np/2+1:Np);
clear S_ac_test_2;
S_ac_test_buling(1:Np/2,:) = S_ac_test_buling_1(1:Np/2,:);
S_ac_test_buling(Up_Sampling_Ratio*Np-(Np-Np/2)+1:Up_Sampling_Ratio*Np,:) =...
    S_ac_test_buling_1(Np/2+1:Np,:);
clear S_ac_test_buling_1;
% figure;imagesc(abs(fftshift(S_ac_test_buling)));title('补零后的二维频谱');
% ========================================================================

S_ac_test_1 = ifft(S_ac_test_buling,[],2);
s_ac_test = ifft(S_ac_test_1,[],1);         % 完成二维升采样
clear S_ac_test_buling;clear S_ac_test_1;

% 作图
% figure;
% imagesc(abs(s_ac_test));
% title('将上述点目标切片做二维升采样');

%%
% （2）读取二维升采样后的峰值点幅相信息
[~,p] = max(abs(s_ac_test));
[~,q] = max(max(abs(s_ac_test)));

% 二维矩阵最大值，所在第几行—— row_max;
row_max = p(q);
% 二维矩阵最大值，所在第几列—— column_max;
column_max = q;
% 矩阵峰值点的幅相信息——data_Peak;
data_Peak = s_ac_test(row_max,column_max);


%%
%
% 计算每个通道下定标点处的信噪比

tmp(1) = var(abs(reshape(s_ac_test(1:1+128-1,1:1+128-1),1,[])));
tmp(2) = var(abs(reshape(s_ac_test(1:1+128-1,end-128+1:end),1,[])));
tmp(3) = var(abs(reshape(s_ac_test(end-128+1:end,1:1+128-1),1,[])));
tmp(4) = var(abs(reshape(s_ac_test(end-128+1:end,end-128+1:end),1,[])));
% max(tmp)        % 选定区域内，噪声方差（无单位）
SNR_tmp = 10*log10(abs(data_Peak)^2/max(tmp));	% 计算信噪比大小，单位 dB

%}

end

