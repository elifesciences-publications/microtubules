%%%% Charlie Jeynes %%%%
%%%% 29/08/2018 %%%%%%%%
%%%% fits fluorescence microtubule data against time with a sigmoidal
%%%% curve, and then performs and ANOVA on the IC50 value (similar to the
%%%% Kd value) of every single data set (in this case 24 data sets) 

clc
clear 
close all

%%
filename = 'For Charlie.xlsx';
T = readtable(filename);

%% create table of fits for each 
Nm =T.Properties.VariableDescriptions; 
x = 1:61;
x =x'; 
param = []; 
for i = 1:24 %number of measurements in the table
    measurement = T{:, i}; 
    param(i,:)  = sigm_fit(x,measurement);
    
        if contains(Nm(i), 'turc', 'IgnoreCase',true) || contains(Nm(i), 'TURC')
            groupNm{i} = 'turc'; 
        elseif contains(Nm(i), 'augmin', 'IgnoreCase',true) 
            groupNm{i} = 'augmin';
        elseif contains(Nm(i), 'WT', 'IgnoreCase',true) 
            groupNm{i} = 'WT';
         elseif contains(Nm(i), 'A', 'IgnoreCase',true) 
            groupNm{i} = 'AT';
        end
    
end
%% perform 1wayANOVA
groupNmT = groupNm'; 
% groupy = repmat(groupNm, 4); 
[~,~, stats] = anova1(param(:, 3), groupNm); 

%% compare means in ANOVA
c = multcompare(stats); 

%% This fits all at the same time
turc = T{:, 1:6}; 
augmin = T{:, 7:12};
WT = T{:, 13:18};
AT = T{:, 19:24};
x = 1:61;
x =x'; 
xx = repmat(x, 1,6); 
figure, 
hold on
[param0,stat0]=sigm_fit(xx,turc);  % param = [min, max, x50, slope]
[param1,stat1]=sigm_fit(xx,augmin); 
[param2,stat2]=sigm_fit(xx,WT); 
[param3,stat3]=sigm_fit(xx,AT); 

%%


function [param,stat]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
    % Optimization of parameters of the sigmoid function
    %
    % Syntax:
    %       [param]=sigm_fit(x,y)       
    %
    %       that is the same that
    %       [param]=sigm_fit(x,y,[],[],[])     % no fixed_params, automatic initial_params
    %
    %       [param]=sigm_fit(x,y,fixed_params)        % automatic initial_params
    %       [param]=sigm_fit(x,y,[],initial_params)   % use it when the estimation is poor
    %       [param]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
    %
    % param = [min, max, x50, slope]
    %
    % if fixed_params=[NaN, NaN , NaN , NaN]        % or fixed_params=[]
    % optimization of "min", "max", "x50" and "slope" (default)
    %
    % if fixed_params=[0, 1 , NaN , NaN]
    % optimization of x50 and slope of a sigmoid of ranging from 0 to 1
    %
    %
    % Additional information in the second output, STAT
    % [param,stat]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
    %
    %
    % Example:
    % %% generate data vectors (x and y)
    % fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)))
    % param=[0 1 5 1];  % "min", "max", "x50", "slope"
    % x=0:0.1:10;
    % y=fsigm(param,x) + 0.1*randn(size(x));
    %
    % %% standard parameter estimation
    % [estimated_params]=sigm_fit(x,y)
    %
    % %% parameter estimation with forced 0.5 fixed min
    % [estimated_params]=sigm_fit(x,y,[0.5 NaN NaN NaN])
    %
    % %% parameter estimation without plotting
    % [estimated_params]=sigm_fit(x,y,[],[],0)
    %
    %
    % Doubts, bugs: rpavao@gmail.com
    % Downloaded from http://www.mathworks.com/matlabcentral/fileexchange/42641-sigmoid-logistic-curve-fit

    % warning off

    x=x(:);
    y=y(:);

    if nargin<=1 %fail
        fprintf('');
        help sigm_fit
        return
    end

    automatic_initial_params=[quantile(y,0.05) quantile(y,0.95) NaN 1];
    if sum(y==quantile(y,0.5))==0
        temp=x(y==quantile(y(2:end),0.5));    
    else
        temp=x(y==quantile(y,0.5));
    end
    automatic_initial_params(3)=temp(1);

    if nargin==2 %simplest valid input
        fixed_params=NaN(1,4);
        initial_params=automatic_initial_params;
        plot_flag=1;    
    end
    if nargin==3
        initial_params=automatic_initial_params;
        plot_flag=1;    
    end
    if nargin==4
        plot_flag=1;    
    end

    if exist('fixed_params','var')
        if isempty(fixed_params)
            fixed_params=NaN(1,4);
        end
    end
    if exist('initial_params','var')
        if isempty(initial_params)
            initial_params=automatic_initial_params;
        end
    end
    if exist('plot_flag','var')
        if isempty(plot_flag)
            plot_flag=1;
        end
    end

    %p(1)=min; p(2)=max-min; p(3)=x50; p(4)=slope como em Y=Bottom + (Top-Bottom)/(1+10^((LogEC50-X)*HillSlope))
    %f = @(p,x) p(1) + (p(2)-p(1)) ./ (1 + 10.^((p(3)-x)*p(4)));

    f_str='f = @(param,xval)';
    free_param_count=0;
    bool_vec=NaN(1,4);
    for i=1:4;
        if isnan(fixed_params(i))
            free_param_count=free_param_count+1;
            f_str=[f_str ' param(' num2str(free_param_count) ')'];
            bool_vec(i)=1;
        else
            f_str=[f_str ' ' num2str(fixed_params(i))];
            bool_vec(i)=0;
        end
        if i==1; f_str=[f_str ' + (']; end
        if i==2;
            if isnan(fixed_params(1))            
                f_str=[f_str '-param(1) )./ (   1 + 10.^( (']; 
            else
                f_str=[f_str '-' num2str(fixed_params(1)) ')./ (1 + 10.^((']; 
            end
        end    
        if i==3; f_str=[f_str ' - xval ) *']; end
        if i==4; f_str=[f_str ' )   );']; end
    end

    eval(f_str)

    [BETA,RESID,J,COVB,MSE] = nlinfit(x,y,f,initial_params(bool_vec==1));
    stat.param=BETA';

    % confidence interval of the parameters
    stat.paramCI = nlparci(BETA,RESID,'Jacobian',J);

    % confidence interval of the estimation
    [stat.ypred,delta] = nlpredci(f,x,BETA,RESID,'Covar',COVB);
    stat.ypredlowerCI = stat.ypred - delta;
    stat.ypredupperCI = stat.ypred + delta;

    % plot(x,y,'ko') % observed data
    % hold on
    % plot(x,ypred,'k','LineWidth',2)
    % plot(x,[lower,upper],'r--','LineWidth',1.5)

    free_param_count=0;
    for i=1:4;
        if isnan(fixed_params(i))
            free_param_count=free_param_count+1;
            param(i)=BETA(free_param_count);
        else
            param(i)=fixed_params(i);
        end    
    end

    if plot_flag==1 
        x_vector=min(x):(max(x)-min(x))/100:max(x);
        plot(x,y,'k.',x_vector,f(param(isnan(fixed_params)),x_vector),'r-')
        xlim([min(x) max(x)])
    end
end

% %% Hill equation
% 
% % Rate = Vmax * S^n / (K^n +S^n)
% 
% hill = @(S, Vmax, K, n) Vmax*S^n + (K^n + S^n); 
% hill_1 = @(S) hill(S,1,0.1,0.5); 
% fplot(hill_1, [0,3])
% %% import data in excel spreadsheet 
% 
% % filename = fullfile(matlabroot,'examples','matlab','myCsvTable.dat');
% 
 
% %% Plot turc
% 
% var = T(:, 1); 
% time = 1:61; 
% figure, 
% plot(time, T{:, 1:6})
% %% example hill equation
% Agonist = time'; 
% atRA = T{:, 1}; 
% % Agonist = [0.1 0.5 1 5 10 19 50 100 114 500 1000 2000];
% % atRA = [0 0 7 15 30 50 58 80 83 87 90 90];
% % MAPPING: Emax = b(1),  EC50 = b(2)
% hill_fit = @(b,x)  b(1).*x./(b(2)+x);
% b0 = [4000; 30];  % b0 = [90; 19];                                  % Initial Parameter Estimates
% B = lsqcurvefit(hill_fit, b0, Agonist, atRA);
% AgVct = linspace(min(Agonist), max(Agonist));   % Plot Finer Resolution
% figure(1)
% plot(Agonist, atRA, 'bp')
% hold on
% plot(AgVct, hill_fit(B,AgVct), '-r')
% hold off
% grid
% xlabel('Agonist')
% ylabel('atRA')
% legend('Data', 'Hill Equation Fit', 'Location','SE')




% ypred = mean(reshape(stat.ypred, 61,6), 6); 
% ypredlowerCI = reshape(stat.ypredlowerCI, 61,6);
% ypredupperCI = mean(reshape(stat.ypredupperCI, 61,6),6);
% figure, 
% hold on
% plot(y)
% plot(ypred, 'o')
% plot(ypredlowerCI, 'x')


% figure, 
% plot(x ,y)
