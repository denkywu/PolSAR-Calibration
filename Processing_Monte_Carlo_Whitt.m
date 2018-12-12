%%
% 蒙特卡洛仿真实验
%
% 考虑噪声 N 的情况下
% 分析不同信噪比下，“检测用三面角反射器”经过极化校准后的极化精度
% 主要包括串扰水平，幅度不平衡和相位不平衡
%
% 更新至：2017.10.28. 17:03


%%
close all;
clear;
clc;


%% 仿真参数设置（根据需要修改）
delta_dB = -500;% -500;-40;-35;-30;-25;-20;     	% 失真矩阵串扰大小，单位 dB

f_dB = 3;% 0;1;2;3;                  % 失真矩阵幅度不平衡，单位 dB


%% 仿真参数设置（不要改！）
f_phase = 0;                    % 失真矩阵相位不平衡，单位 度（°）
SNR_dB = 45 : -1 : 25;          % 极化定标模型中信噪比大小，单位 dB

Num_Monte_Carlo = 10000;      % 蒙特卡洛实验仿真次数（某一SNR下）


%% 仿真
h = waitbar(0,'Please wait...');
Flag = 0; % Flag == 0，表示利用 22.5°二面角反射器进行Whitt算法求解（勿动！）

for qq = 1 : length(SNR_dB)
    for pp = 1 : Num_Monte_Carlo
        % 生成数据
        [T,R,S1,S2,S3,S_test,M1,M2,M3,M_test,M_test2] = Gen_Data_Whitt(delta_dB,f_dB,f_phase,SNR_dB(qq));

        % 用Whitt算法求解失真矩阵
        [T_solve,R_solve,A_solve] = My_Whitt(S1,M1,S2,M2,S3,M3,M_test,Flag);
        
        % 对"检测用三面角反射器"进行极化处理
        M_After_PC = M_Polarization(reshape(M_test2.',1,[]),T_solve,R_solve,A_solve);

        % 将仿真结果进行汇总,便于后续分析
        %   1）zhenzhi_T_R_A 表示每次仿真中失真矩阵的真值
        zhenzhi_T_R_A(pp,1:4) = reshape(T.',1,[]);          % T
        zhenzhi_T_R_A(pp,4+1:4+4) = reshape(R.',1,[]);      % R
        zhenzhi_T_R_A(pp,4+4+1) = 1;                        % A，始终设置为1
        %   2）Resolve_T_R_A 表示用Whitt算法求解得到的失真矩阵，求解值
        Resolve_T_R_A(pp,1:4) = reshape(T_solve.',1,[]);    % T_solve
        Resolve_T_R_A(pp,4+1:4+4) = reshape(R_solve.',1,[]);% R_solve
        Resolve_T_R_A(pp,4+4+1) = A_solve;                  % A_solve
        %   3）Result_M_test 表示“检测用三面角反射器”的极化校准结果
        Result_M_test(pp,1:4) = M_After_PC;

        waitbar(((qq-1)*Num_Monte_Carlo+pp)/(length(SNR_dB)*Num_Monte_Carlo));
    end

    tmp = Result_M_test(:,1);
    guiyi_Result_M_test = Result_M_test./(tmp*ones(1,4));   % 对 Result_M_test 归一化
    clear tmp;
    
    %% 统计该SNR下的“最差指标”
    dB_Result_M_test = 20*log10(abs(guiyi_Result_M_test));  % 取幅度（dB）
    zuicha_HV_VH(qq,1) = max(max(dB_Result_M_test(:,2:3))); % 串扰（dB），HV和VH通道一起分析
    zuicha_f_abs(qq,1) = max(abs(dB_Result_M_test(:,4)));   % 幅度不平衡（dB），这是在正负的基础上取绝对值得到的最大。也就是说幅度不平衡在 [-XX，+XX]之间。
    angle_guiyi_Result_M_test = angle(guiyi_Result_M_test)/pi*180;% 取相位（°） 
    zuicha_f_phase(qq,1) = max(abs(angle_guiyi_Result_M_test(:,4)));% 相位不平衡（dB），这是在正负的基础上取绝对值得到的最大。类似幅度不平衡。
    clear dB_Result_M_test;
    
    %% 统计该SNR下的“（μ+σ）指标”
    abs_guiyi_Result_M_test = abs(guiyi_Result_M_test);% 取幅度（非dB）
    temp = abs_guiyi_Result_M_test(:,2:3);
    temp = reshape(temp,[],1);
    miusigma_HV_VH(qq,1) = 20*log10( mean(temp) + 1 * std(temp) );% 串扰不符合正太分布，用（μ+σ）；单位dB
    clear temp;
    miusigma_f_abs(qq,1) = 20*log10( mean(abs_guiyi_Result_M_test(:,4)) + 2 * std(abs_guiyi_Result_M_test(:,4)) );% 幅度近似正太分布，用（μ+2σ），95%置信区间；单位dB
    miusigma_f_phase(qq,1) = ( mean(angle_guiyi_Result_M_test(:,4)) + 2 * std(angle_guiyi_Result_M_test(:,4)) );% 相位不平衡近似正太分布，用（μ+2σ），95%置信区间；单位度
   
    clear angle_Result_M_test;
    clear abs_guiyi_Result_M_test;
    clear guiyi_Result_M_test;
    
end
close(h);
clc;


%%
% 当 SNR_dB 选择为一组不同强度的信噪比时
% 分析“检测用三面角反射器”极化处理结果随信噪比的变化
%% 极化精度 1——最差指标
%
figure;plot(SNR_dB,zuicha_HV_VH);grid on;
title('最差指标，串扰大小随信噪比的变化');
xlabel('信噪比，单位 dB');ylabel('串扰大小，单位 dB');legend('无串扰,f=3dB');%legend('delta=-20dB,f=0dB');
figure;plot(SNR_dB,zuicha_f_abs);grid on;
title('最差指标，幅度不平衡随信噪比的变化');
xlabel('信噪比，单位 dB');ylabel('幅度不平衡，单位 dB');legend('无串扰,f=3dB');
figure;plot(SNR_dB,zuicha_f_phase);grid on;
title('最差指标，相位不平衡随信噪比的变化');
xlabel('信噪比，单位 dB');ylabel('相位不平衡，单位 度');legend('无串扰,f=3dB');
%}
%% 极化精度 2——（μ+σ）指标
%
figure;plot(SNR_dB,miusigma_HV_VH);grid on;
title('（μ+σ）指标，串扰大小随信噪比的变化');
xlabel('信噪比，单位 dB');ylabel('串扰大小，单位 dB');legend('无串扰,f=3dB');
figure;plot(SNR_dB,miusigma_f_abs);grid on;
title('（μ+2σ）指标，幅度不平衡随信噪比的变化');
xlabel('信噪比，单位 dB');ylabel('幅度不平衡，单位 dB');legend('无串扰,f=3dB');
figure;plot(SNR_dB,miusigma_f_phase);grid on;
title('（μ+2σ）指标，相位不平衡随信噪比的变化');
xlabel('信噪比，单位 dB');ylabel('相位不平衡，单位 度');legend('无串扰,f=3dB');
%}


%% 某一SNR下的分析
%{
% “检测用三面角反射器”极化处理结果
tmp = Result_M_test(:,1);
guiyi_Result_M_test = Result_M_test./(tmp*ones(1,4));   % 对 Result_M_test 归一化
clear tmp;

dB_Result_M_test = 20*log10(abs(guiyi_Result_M_test));  % 幅度，单位 dB
figure;plot(dB_Result_M_test(:,2),'.');
title('HV 串扰大小');xlabel('仿真次数');ylabel('大小，单位 dB')
hold on;plot(max(dB_Result_M_test(:,2)).*ones(Num_Monte_Carlo,1),'-r','LineWidth',2);hold off;
figure;plot(dB_Result_M_test(:,3),'.');
title('VH 串扰大小');xlabel('仿真次数');ylabel('大小，单位 dB')
hold on;plot(max(dB_Result_M_test(:,3)).*ones(Num_Monte_Carlo,1),'-r','LineWidth',2);hold off;
figure;plot(dB_Result_M_test(:,4),'.');
title('幅度不平衡');xlabel('仿真次数');ylabel('大小，单位 dB')
hold on;
plot(max(dB_Result_M_test(:,4)).*ones(Num_Monte_Carlo,1),'-r','LineWidth',2);
plot(min(dB_Result_M_test(:,4)).*ones(Num_Monte_Carlo,1),'-r','LineWidth',2);
hold off;

angle_guiyi_Result_M_test = angle(guiyi_Result_M_test)/pi*180;% 相位，单位度（°） 
figure;plot(angle_guiyi_Result_M_test(:,4),'.');
title('相位不平衡');xlabel('仿真次数');ylabel('相位，单位 度（°）');
hold on;
plot(max(angle_guiyi_Result_M_test(:,4)).*ones(Num_Monte_Carlo,1),'-r','LineWidth',2);
plot(min(angle_guiyi_Result_M_test(:,4)).*ones(Num_Monte_Carlo,1),'-r','LineWidth',2);
hold off;

%%
% 矩阵 T 的真值和求解值
figure;
suptitle('失真矩阵 T 的求解值和真值对比')

subplot(2,2,1);
plot(20*log10(abs(Resolve_T_R_A(:,2)./Resolve_T_R_A(:,1))))
title('T 的 HV 串扰大小');xlabel('仿真次数');ylabel('大小，单位 dB')
hold on;plot(20*log10(abs(zhenzhi_T_R_A(:,2))),'r-','LineWidth',2);hold off;
subplot(2,2,2);
plot(20*log10(abs(Resolve_T_R_A(:,3)./Resolve_T_R_A(:,1))))
title('T 的 VH 串扰大小');xlabel('仿真次数');ylabel('大小，单位 dB')
hold on;plot(20*log10(abs(zhenzhi_T_R_A(:,3))),'r-','LineWidth',2);hold off;
subplot(2,2,3);
plot(20*log10(abs(Resolve_T_R_A(:,4)./Resolve_T_R_A(:,1))))
title('T 的幅度不平衡');xlabel('仿真次数');ylabel('大小，单位 dB')
hold on;plot(20*log10(abs(zhenzhi_T_R_A(:,4))),'r-','LineWidth',2);hold off;
subplot(2,2,4);
plot((angle(Resolve_T_R_A(:,4)./Resolve_T_R_A(:,1)))/pi*180)
title('T 的相位不平衡');xlabel('仿真次数');ylabel('相位，单位 度（°）');
hold on;plot((angle(zhenzhi_T_R_A(:,4)))/pi*180,'r-','LineWidth',2);hold off;

%%
% 矩阵 R 的真值和求解值
figure;
suptitle('失真矩阵 R 的求解值和真值对比')

subplot(2,2,1);
plot(20*log10(abs(Resolve_T_R_A(:,4+2)./Resolve_T_R_A(:,4+1))))
title('R 的 HV 串扰大小');xlabel('仿真次数');ylabel('大小，单位 dB')
hold on;plot(20*log10(abs(zhenzhi_T_R_A(:,4+2))),'r-','LineWidth',2);hold off;
subplot(2,2,2);
plot(20*log10(abs(Resolve_T_R_A(:,4+3)./Resolve_T_R_A(:,4+1))))
title('R 的 VH 串扰大小');xlabel('仿真次数');ylabel('大小，单位 dB')
hold on;plot(20*log10(abs(zhenzhi_T_R_A(:,4+3))),'r-','LineWidth',2);hold off;
subplot(2,2,3);
plot(20*log10(abs(Resolve_T_R_A(:,4+4)./Resolve_T_R_A(:,4+1))))
title('R 的幅度不平衡');xlabel('仿真次数');ylabel('大小，单位 dB')
hold on;plot(20*log10(abs(zhenzhi_T_R_A(:,4+4))),'r-','LineWidth',2);hold off;
subplot(2,2,4);
plot((angle(Resolve_T_R_A(:,4+4)./Resolve_T_R_A(:,4+1)))/pi*180)
title('R 的相位不平衡');xlabel('仿真次数');ylabel('相位，单位 度（°）');
hold on;plot((angle(zhenzhi_T_R_A(:,4+4)))/pi*180,'r-','LineWidth',2);hold off;

%%
% A 的真值和求解值
figure;
plot(Resolve_T_R_A(:,4+4+1));
title('绝对幅度因子 A 的大小');xlabel('仿真次数');
hold on;plot(zhenzhi_T_R_A(:,4+4+1),'r-','LineWidth',2);hold off;


%}
