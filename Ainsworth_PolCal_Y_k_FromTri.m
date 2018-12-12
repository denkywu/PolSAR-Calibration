function [Y,k] = Ainsworth_PolCal_Y_k_FromTri(u,v,w,z,alpha,O)
% 基于三面角反射器数据，求解 Ainsworth 参数 Y 和 k。
% 注1： Ainsworth 原始文献中的下标是“前收 后发”；
%       而室里的下标都是“前发 后收”；
% 注2：该函数的输入是室里规则，因此在实现时需要首先对观测矩阵O的顺序做调整；
%
% 输入：
%   1）已定标参数：u,v,w,z,alpha
%   2）三面角反射器的观测矩阵 O = [ Ohh,Ohv,Ovh,Ovv ];% 这是【室里规则】“前发 后收”！
%
% 输出：
%   定标参数 Y 和 k ;
%
% 本程序更新至：2017.12.25. 15:51


%%
% 注：由于我的输入是按照室里规则，“前发 后收”表示的；
%     而 Ainsworth 定标模型所基于的下标是“前收 后发”；
%     因此需要对观测矩阵O的顺序做调整；
O2 = [ O(1),O(3),O(2),O(4) ];% 重新排列顺序：这样的下标才对应Ainsworth算法。
clear O;
O = O2; % 行向量，O = [ Ohh,Ohv,Ovh,Ovv ]，该下标是国际惯例（前收 后发）。
clear O2;


%%
A = [   alpha^2,        v,          w*alpha^2,      v*w;
        z*alpha^2,      1,          w*z*alpha^2,    w;
        u*alpha^2,      u*v,        alpha^2,        v;
        u*z*alpha^2,    u,          z*alpha^2,      1   ];

S_O = A\(O.');          % 也即 S_O = inv(A)*(O.');
                       	% 这里根据 MATLAB 提示将 INV(A)*b 修改为 A\b;  
% S_O 理论上等于下式
%       S_O = Y.*[  k^2*Shh;
%                   k*Shv;
%                   k*Svh;
%                   Svv     ];% 下标国际惯例


%%
% 下面利用 S_O 求解 Y 和 k 的大小
%   1）求解 Y
Y = S_O(4);

%   2）求解 k
arg_k = angle( S_O(1)/Y )/2;                    % k 的相位（会有180°的模糊）
abs_k = ( abs(S_O(1))^2/(abs(Y)^2) )^(1/4);  	% k 的幅度
k = abs_k*exp(1j*arg_k);


%%
% 至此，Y 和 k 的数值已求解得到。


end