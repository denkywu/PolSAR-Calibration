function M_After_PC = M_Polarization(M,T_solve,R_solve,A_solve)
% 用求解得到的失真矩阵对测试定标点进行定标处理,由观测矩阵反解求得测试定标点的散射矩阵
% 具体方法：
%   由 M = A*RST 推得：S = (1/A)*R^(-1)*M*T^(-1);
%   其中 M 表示测试定标点的观测矩阵，也采用峰值法提取;
%   求解得到 S（M_After_PC）,与理论散射矩阵比较即可进行指标测试
%
% 输入变量：
%   1）M             观测矩阵（测试用的普通三面角反射器的，非极化定标器）
%   2）T_solve       发射失真矩阵
%   3）R_solve       接收失真矩阵
%   4）A_solve       绝对幅度因子
%
% 输出变量：
%   1）M_After_PC    定标后的结果（求解得到的散射矩阵 S）
%
% 该程序截止至：2016.05.24. 15:53

%%
% 用测试用三面角反射器进行指标测试
[num_M,~] = size(M);            % M中每一行表示不同的定标器，总共 num_M 个;
                                % 每一行的4个值依次代表 HH,HV,VH,VV 通道;
                                
M_After_PC = zeros(num_M,4);	% 用来存放定标后的结果：
                                % 每一行表示一个定标器的处理结果;
                                % 每一行的4个值依次代表 HH,HV,VH,VV 通道;
for pp = 1:num_M
    M_test = M(pp,:);
    M_test = reshape(M_test,[2,2]).';
    S_test = (1/A_solve)*R_solve^(-1)*M_test*T_solve^(-1);
    M_After_PC(pp,:) = reshape(S_test.',1,[]);
    
    clear M_test;clear S_test;
end

% disp('-------------------------------------------------------------')
% disp('测试指标结果，计算完成');
% disp('-------------------------------------------------------------')

end