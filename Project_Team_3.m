%% Group project team 3: compute the VAR for Operational losses
clear 
close all
clc
%% Data analysis
% Graphic of data
load LossRisk.mat;
losses=oprisk;
figure(1);
histogram(losses,100);
ylim([0 250]);
xlim([0 2.3*10^6])
xlabel('Losses','FontSize',12)
ylabel('Frequency','FontSize',12)
title('Operational losses','FontSize',14)

% Descriptive Statistics
Stats.Mean_Losses=mean(losses);
Stats.Std_Losses=std(losses);
Stats.Median_Losses=median(losses);
Stats.Skew_Losses=skewness(losses);
Stats.Kurt_Losses=kurtosis(losses);
Stats.Max_Losses=max(losses);
Stats.Min_Losses=min(losses);
disp(Stats)

% Var on empirical data
alpha=[0.05 0.01 0.005 0.0001];
for i=1:4
    VarEmp(i)=quantile(losses,1-alpha(i));
end

%% Lognormal distribution approach
lognorm=fitdist(losses,"Lognormal");
logv=logncdf(losses,lognorm.ParameterValues(1,1),lognorm.ParameterValues(1,2));
figure(2)
histfit(losses,100,"lognormal");
ylim([0 500])
xlim([0 2.3*10^6])
xlabel('Losses','FontSize',12)
ylabel('Frequency','FontSize',12)
title('Fitted lognormal distribution on operational losses','FontSize',14)

%% test fit distributions
testLogn=makedist("Lognormal","mu",lognorm.ParameterValues(1,1),"sigma",lognorm.ParameterValues(1,2));
[H,pvalueLog,ksstatlog]=kstest(losses,"Alpha",0.05,"CDF",testLogn);
tabLogpar=table(lognorm.ParameterValues(1,1),lognorm.ParameterValues(1,2),pvalueLog,VariableNames=["Mu","Sigma","KStest_pvalue"]);
disp(tabLogpar)
%% Lognormal Var
alpha=[0.05 0.01 0.005 0.0001];
for i=1:4
VarLognorm(i)=logninv(1-alpha(i),lognorm.ParameterValues(1,1),lognorm.ParameterValues(1,2));
end

%% EVT for operational risk
% selection of the threshold
u_startingvalue=quantile(oprisk,0.90);
u_finalvalue=quantile(oprisk,0.98);
thresh_seq = 100000:10000:318347; % choose a threshold, compute the MEF for various level of u
losses=oprisk;
for i=1:length(thresh_seq)
    excesses = losses(losses > thresh_seq(i)) - thresh_seq(i); % set the excesses for each u
    Nuv(i) = length(excesses);
    mef(i) = sum(excesses)/Nuv(i); % emprical mef
end
figure(3)
plot(thresh_seq,mef)
xlabel('Treshold value','FontSize',12)
ylabel('MEF value','FontSize',12)
title('MEF function','FontSize',14)
% we set the threshold
thresh=220000; % at this value the mef function tends to be linear (proxy)
% We calculate the excesses for the level of thresh selected
excesses = losses(losses > thresh) - thresh; % set the excesses for u
exceedances = losses(losses > thresh); % optional( needed for the graph)
% estimation of the parameter of the GPD with the POT method
n=length(oprisk);
Nu = length(excesses)/n;
[parmhat,parmci] = gpfit(excesses); % shape, scale confirms that is a GPD
% Graphic of the Pareto distribution
exc=thresh:1000:max(losses);
gpdpdf = gppdf(exc,parmhat(1),parmhat(2),thresh);
figure(4)
histogram(exceedances,100,'normalization','pdf')
hold on
d = plot(exc,gpdpdf,'k','LineWidth',2);
xlabel('Exceedances','FontSize',12)
ylabel('Frequency','FontSize',12)
title('Fitted Generalized Pareto Distribution on operational losses','FontSize',14)
% Test KS
pd=fitdist(excesses,'GeneralizedPareto');
[h,pvalueGpd,ksstatGpd]=kstest(excesses,'CDF',pd);
% We put in a table the parameters of the GPD distribution
tabGPDpar=table(parmhat(1),parmhat(2),pvalueGpd,VariableNames=["Shape","Scale","KStest_pvalue"])
% p-value=0.3235 it is higher than 5%, so we don't reject the null hypotesis
% estimate the VAR
alpha=[0.05; 0.01; 0.005; 0.0001];
Nexc=length(excesses);
for i=1:length(alpha)
    VarEvt(i)=thresh + (parmhat(2)/parmhat(1))*((n*alpha(i)/Nexc)^(-parmhat(1))-1);
end
% Compare the Var of the two model and put them in table
tabVar=table(alpha,VarEmp',VarLognorm',VarEvt',VariableNames=["Alpha","VarEmp", "VarLognorm","VarEvt"]);
disp(tabVar)
% estimates tail probability for different level of the threshold, in order
% to assess which model is more suitable for LFHI events. Here we compare
% the probability with a GPD and the empirical probability from the data
% (need to add the lognormal probability)
mu=lognorm.ParameterValues(1,1);
sigma=lognorm.ParameterValues(1,2);
x=[200000; 250000; 300000; 350000; 400000; 450000; 577000; 2000000; 10000000];
for i=1:length(x)
    observed(i)=length(losses(losses>x(i)))/n;
    p_LogN(i)=1-logncdf(x(i),mu,sigma);
    p_EVT(i)=Nu*(1+parmhat(1)*(x(i)-thresh)/parmhat(2))^(-1/parmhat(1));
end
% we put the variables in a table
Threshold={'200000'; '250000'; '300000'; '350000'; '400000'; '450000'; '577000'; '2000000'; '10000000'};
observed=observed';
p_LogN=p_LogN';
p_EVT=p_EVT';
tabProb=table(Threshold,observed,p_LogN,p_EVT);
disp(tabProb)
% Graphical comparison between the LogNormal and the GPD, which model is
% the best?
figure(5)
histogram(exceedances,100,'normalization','pdf')
hold on
l=plot(exc,TruncLogN(exc,mu,sigma,[thresh Inf]),'r','LineWidth',2);
e = plot(exc,gpdpdf,'k','LineWidth',2);
xlabel('Exceedances','FontSize',12)
ylabel('Frequency','FontSize',12)
title('Operational Losses fitted by a LogN and a GPD densities','FontSize',14)
legend([l e],'LogNormal','GPD','Location','northeast','Fontsize',12)

%% Case with the trashold at 20000
% we set the threshold arbitrary at 20000 because with that value the EVT
% approach give largere Var at 99.99%
thresh1=20000;
% We calculate the excesses for the level of thresh selected
excesses = losses(losses > thresh1) - thresh1; % set the excesses for u
exceedances = losses(losses > thresh1); % optional( needed for the graph)
% estimation of the parameter of the GPD with the POT method
n=length(losses);
Nu = length(excesses)/n;
[parmhat,parmci] = gpfit(excesses); % shape, scale confirms that is a GPD
% Graphic of the Pareto distribution
exc=thresh:1000:max(losses);
gpdpdf1 = gppdf(exc,parmhat(1),parmhat(2),thresh1);
% Test KS
pd=fitdist(excesses,'GeneralizedPareto');
[h,pvalueGpd1,ksstatGpd1]=kstest(excesses,'CDF',pd);
% We put in a table the parameters of the new GPD distribution
tabGPDpar=table(parmhat(1),parmhat(2),pvalueGpd1,VariableNames=["Shape","Scale","KStest_pvalue"])
% p-value=0.1921 is higher than 5%, so we don't reject the null hypotesis
% estimate the VAR
alpha=[0.05; 0.01; 0.005; 0.0001];
Nexc=length(excesses);
for i=1:length(alpha)
    VarEvt1(i)=thresh1 + (parmhat(2)/parmhat(1))*((n*alpha(i)/Nexc)^(-parmhat(1))-1);
end
% Compare the Var of the two model and put them in table
tabVar1=table(alpha,VarLognorm',VarEvt1',VariableNames=["Alpha", "VarLognorm","VarEvt"]);
disp(tabVar1)


