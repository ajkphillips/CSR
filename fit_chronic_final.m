%% Code for fitting chronic sleep restriction model to Van Dongen data.

load vandongenbu % Load data from Van Dongen et al. chronic sleep restriction study

pars0 = [220, 8.4, 0.9]; %,300]; %, 0.0198];
%pars0 = [436.16, 10.039, 0.94015]; %, 19.44];
times = [dtotal(:,1);dp4(:,1);dp6(:,1);dp8(:,1)]; % Times at which PVT was measured
PVT_exp = [dtotal(:,2);dp4(:,2);dp6(:,2);dp8(:,2)]; % PVT lapses from experiment

%beta0 = [3.5,100,0.45,69,656,63,35.2,7.9,20,0.005]; % initial guess at parameter values
%beta0 = 500;
beta0 = [1, 100, 0.85, 581.6, 6.5, 3.2, 8.25]; % default initial guess for fit in paper
%beta0 = [3.5, 100, 0.5, 600, 6, 3, 9]; 
%beta0 = [1,100,0.85,579,2.8,1.3,8.2]; % fit from going down to [0.5,50]
%beta0 = [1,100,0.90,565,0.23,-0.09,20.4]; % fit from going down to [0.1,10]
%beta0 = [0.09,9,0.91,563,0.19,0.07,32]; % fit from going down to [0.09,9]
%beta0 = [log10(0.05),log10(5),0.91,563,0.11,0.04,56]; % fit from going down to [0.05,4]
%beta0 = [0.04,5,0.90005,542.54,0.076695,0.029001,56.391];
%beta0 = [0.014991,5.008,0.90067,542.77,0.028617,0.010826,56.388]; % iteration limit exceeded, decreasing only Kd1
%beta0 = [0.013227,5.0104,0.90096,542.95,0.025191,0.0095341,56.387]; % next limit exceeded, allowing 0.01, 1
%beta0 = [0.012865,5.0106,0.90098,542.96,0.024511,0.0092707,56.387];
%beta0 = [0.012529,5.0109,0.901,542.97,0.023868,0.0090275,56.387];
%beta0 = [0.012213,5.0112,0.90102,542.98,0.023263,0.0087986,56.387];
%beta0 = [0.01,5.0112,0.90105,543,0.023,0.0085,56.387]; % fresh guess to try to go lower
%beta0 = [0.011286,5.0119,0.90108,543.01,0.021488,0.0081269,56.387];
%beta0 = [0.011067,5.0123,0.90111,543.02,0.021068,0.0079681,56.387];
%beta0 = [0.011034,5.0122,0.90111,543.02,0.021004,0.0079438.56.387] % ACTUAL BEST FIT!!
    % guesses for Kd1, Kd2, musleep/Atot, Dmid, Ds, a, phi, 

size(chronic_model_final(beta0,times))

PVT_model = @(par,t) chronic_model_final(par,t);

beta = nlinfit(times,PVT_exp,PVT_model,beta0);

PVTmod = chronic_model_final(beta,times);

figure(24)
plot(PVT_exp)
hold on
plot(PVTmod,'r')
hold off
ylabel('Lapses')
legend('Exp','Mod')

SStot = sum((PVT_exp-mean(PVT_exp)).^2);
SSres = sum((PVT_exp-PVTmod).^2);

R2 = 1 - SSres/SStot;

p = length(beta0);
n = length(PVT_exp);
R2a = R2 - (1-R2)*p/(n-p-1);

disp(['Adjusted R^2 = ',num2str(R2a)])

