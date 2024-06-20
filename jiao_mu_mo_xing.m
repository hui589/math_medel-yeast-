time=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36];
popsize=[0.200,0.330,0.500,1.10,1.40,3.10,3.50,9.00,10.0,25.4,27.0,55.0,76.0,115,160,162,190,193,190,209,190,210,200,215,220,200,180,213,210,210,220,213,200,211,200,208,230];
scatter(time,popsize);
xlabel('时间 (小时)');
ylabel('酵母数量');
%------------------------------------------------------%
ln_popsize=log(popsize);	            %取对数
p = polyfit(time,ln_popsize,1);	        %一次线性拟合得到斜率和截距
r=p(1);			                        %斜率
ln_N0=p(2);		                        %截距
exp_model=@(t)exp(ln_N0 + r * t);   	%增长模型
hold on;
fplot(exp_model,[0 16]);		        %画出图形，时间范围0~16
%------------------------------------------------------%
logistic_ode=@(t,N,r,K)r*N*(1-N/K);	    %定义 Logistic 增长模型的微分方程
r_guess = 0.1;				            %初始猜测的增长率
K_guess = 200;		    	            %初始猜测的环境承载能力
params_guess = [r_guess, K_guess];	    %将猜测的参数放在一个向量中
yeast_count = popsize;%将popsize复制给yeast_count
TN = @(params, time, yeast_count, logistic_ode) ...
    ode45(@(t, N) logistic_ode(t, N, params(1), params(2)), ...
    [0, max(time)], yeast_count(1));                            %计算微分方程
compute_error = @(params, time, yeast_count, logistic_ode) ...
    sum(((interp1(TN(params, time, yeast_count, logistic_ode).x,...
    TN(params, time, yeast_count, logistic_ode).y, time)) - yeast_count).^2);% 计算误差
error_function = @(params) compute_error(params, time, yeast_count, logistic_ode);% 定义误差函数

optimal_params = fminsearch(error_function, params_guess);% 使用 fminsearch 优化参数
r_opt = optimal_params(1); % 最优增长率
K_opt = optimal_params(2); % 最优环境承载能力

[T, N] = ode45(@(t, N) logistic_ode(t, N, r_opt, K_opt), [0, max(time)], yeast_count(1));% 使用优化后的参数求解微分方程
plot(T, N, 'b', 'LineWidth', 2);% 绘制 Logistic 增长模型曲线
legend('数据点', '指数增长模型','Logistic 增长模型');