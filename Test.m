%% Recreating Results from Previous Publications
% [1] https://www.sciencedirect.com/science/article/pii/S1359836813005337
% [2] https://www.sciencedirect.com/science/article/pii/S1110016822008420#b0170


set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

figure
hold on
E1 = 19; E2 = 2; G12 = 1; v12 = 0.3; v21 =  v12*E2/E1;
thickness = 0.0052;
materialList = [E1 E2 G12 v12 thickness];
orientations = linspace(0,90,91);

% my code
E_modulus = zeros(1,length(orientations));
for i = 1:length(orientations) 
    obj = laminateCLT_700(materialList,orientations(i));
    obj = obj.properties;
    E_modulus(i) = obj.myE1;
end
plot(orientations,E_modulus/E1,'r') 

% results from equation in reference [1]
E_1 = E1.*cosd(orientations).^4 + E2.*sind(orientations).^4 + G12.*sind(2.*orientations).^2;
plot(orientations,E_1/E1,'b')

% results from equation in reference [2]
n = 1/E1.*(cosd(orientations).^4 - v12.*sind(orientations).^2.*cosd(orientations).^2)...
    + 1/E2.*(sind(orientations).^4 - v21.*sind(orientations).^2.*cosd(orientations).^2)...
    + 1/G12.*(sind(orientations).^2.*cosd(orientations).^2);
Ex = 1./n;
plot(orientations,Ex/E1,'k--')

legend('Classical Laminate Theory Code', 'Reference [1]', 'Reference [2]')
xlabel('Orientation, $\theta$ (deg)'); ylabel('Normalized $E_1$')
set(gcf,'color','white')


%% Generating Monte Carlo Data of 1-5 layers
% using CLT code

layers = 5;
sample = 1e6;
variability = [1 2 3];

E_normalized = cell(1,3);
E_modulus = zeros(layers,sample);
z = zeros(layers,sample);
theta = zeros(layers,sample);
orientations = cell(1,layers);

tic
for k = 1:3
    for t = 1:layers
        z(t,:) = randn(1,sample);
        theta(t,:) = 0 + variability(k)*z(t,:);
        for q = 1:t
            orientations{t}(q,:) = theta(q,:);
        end
        for i = 1:sample 
        obj = laminateCLT_700(materialList,orientations{t}(:,i));
        obj = obj.properties;
        E_modulus(t,i) = obj.myE1/E1;
        end
    end
    E_normalized{k} = E_modulus;
end
computationTime = toc;

figure
histogram(orientations{1},'Normalization','pdf')

%% Monte Carlo Data Characterization

avg = zeros(layers,3);
standard_dev = zeros(layers,3);
median = zeros(layers,3);
int_95 = zeros(layers,3);
int_99 = zeros(layers,3);

failure_condition = 0.95;
pf = zeros(layers,3);

for k = 1:3
    for t = 1:layers
        E_normalized{k}(t,:) = sort(E_normalized{k}(t,:));
        data = E_normalized{k}(t,:);
        
        avg(t,k) = mean(data);
        standard_dev(t,k) = std(data);
        median(t,k) = data(round(length(data)/2));
        int_95(t,k) = data(round(length(data)/20));
        int_99(t,k) = data(round(length(data)/100));
        pf(t,k) = sum(data < failure_condition)/numel(data)*100;
    end
end

T_mean = array2table(avg);
T_mean.Properties.VariableNames = {'std_1', 'std_2', 'std_3'};
T_mean.Properties.RowNames = {'1 Plies','2 Plies','3 Plies','4 Plies' ,'5 Plies'};
fprintf('Mean Values\n')
disp(T_mean)
figure
for i = 1:3
    plot(1:1:5,avg(:,i))
    hold on
end
title('Mean Values')
xlabel('Ply Count'); ylabel('Normalized $E_1$')
legend('$\sigma_1$','$\sigma_2$','$\sigma_3$')
set(gcf,'color','white')

T_std = array2table(standard_dev);
T_std.Properties.VariableNames = {'std_1', 'std_2', 'std_3'};
T_std.Properties.RowNames = {'1 Plies','2 Plies','3 Plies','4 Plies' ,'5 Plies'};
fprintf('Standard Deviation Values\n')
disp(T_std)
figure
for i = 1:3
    plot(1:1:5,standard_dev(:,i))
    hold on
end
title('Standard Deviation')
xlabel('Ply Count'); ylabel('Normalized $E_1$')
legend('$\sigma_1$','$\sigma_2$','$\sigma_3$')
set(gcf,'color','white')

T_med = array2table(median);
T_med.Properties.VariableNames = {'std_1', 'std_2', 'std_3'};
T_med.Properties.RowNames = {'1 Plies','2 Plies','3 Plies','4 Plies' ,'5 Plies'};
fprintf('Median Values\n')
disp(T_med)
figure
for i = 1:3
    plot(1:1:5,median(:,i))
    hold on
end
title('Median Values')
xlabel('Ply Count'); ylabel('Normalized $E_1$')
legend('$\sigma_1$','$\sigma_2$','$\sigma_3$')
set(gcf,'color','white')

T_95 = array2table(int_95);
T_95.Properties.VariableNames = {'std_1', 'std_2', 'std_3'};
T_95.Properties.RowNames = {'1 Plies','2 Plies','3 Plies','4 Plies' ,'5 Plies'};
fprintf('95%% Confidence Interval Lower Bound\n')
disp(T_95)

T_99 = array2table(int_99);
T_99.Properties.VariableNames = {'std_1', 'std_2', 'std_3'};
T_99.Properties.RowNames = {'1 Plies','2 Plies','3 Plies','4 Plies' ,'5 Plies'};
fprintf("99%% Confidence Interval Lower Bound\n")
disp(T_99)

T_pf = array2table(pf);
T_pf.Properties.VariableNames = {'std_1', 'std_2', 'std_3'};
T_pf.Properties.RowNames = {'1 Plies','2 Plies','3 Plies','4 Plies' ,'5 Plies'};
fprintf('Probability of Failure (<95%% E_1)\n')
disp(T_pf)
figure
for i = 1:3
    plot(1:1:5,pf(:,i))
    hold on
end
title('Probability of Failure')
xlabel('Ply Count'); ylabel('Normalized $E_1$')
legend('$\sigma_1$','$\sigma_2$','$\sigma_3$')
set(gcf,'color','white')

%% Histograms of Monte Carlo Data
% I can use histogram function to make probability, pdf, and cdf graphs.
% These are all probability, because it's the easiest for me to read.

bin_count = 5*10*round((1 + 3.3*log10(sample))/10);
edges = linspace(0,1,bin_count+1);

for k = 1:layers 
    figure
    for t = 1:3
        histogram(E_normalized{t}(k,:),edges,'Normalization','probability')
        hold on
    end
    ylabel('Probability'); xlabel('Normalized $E_1$')
    legend('$\sigma_1$','$\sigma_2$','$\sigma_3$')
    set(gcf,'color','white')
    xlim([0.6 1])
    title([num2str(k) ' Layers'])
end

%% Matching to a Distribution
% plotting the most extreme distribution, data = E_normalized{3}(1,:)
% these plots are matching pdf, not probability

data = E_normalized{3}(1,:);

omega = 1;                      % upper bound 
y = omega - data;               % transform to Wiebell (from reverse Weibell)
[paramEsts, ~] = wblfit(y);     % MATLAB estimate of Wiebell parameters
         
lambda = paramEsts(1);          % scale parameter 
k = paramEsts(2);               % shape parameter
x = 0:0.001:0.997;              % not going to one avoids major spike

z = (omega - x) / lambda;

pdf = (k / lambda) .* z.^(k-1) .* exp(-(z.^k));

figure
plot(x, pdf, 'k')
xlabel('Normalized $E_1$')
ylabel('Probability Density')
title('Weibull Distribution of One Layer, $\sigma_1$')
grid on
hold on
histogram(data,edges,'Normalization','pdf','FaceAlpha',0.4)
hold off
legend('Weibull PDF','Laminate $E_3$ PDF')
set(gcf,'color','white')
xlim([0.5 1])

x = 0:0.001:1;
cdf = exp(-((omega - x)./lambda).^k);

figure
plot(x, cdf,'k')
xlabel('Normalized E_1')
ylabel('Cumulative Density')
title('Weibull Distribution (bounded at \omega)')
grid on
hold on
histogram(data,edges,'Normalization','cdf','FaceAlpha',0.4)
hold off
legend('Weibull CDF','Laminate E_1 CDF')
set(gcf,'color','white')

% chi-square test


% KS test


%% FORM

% demonstrate how FORM would be difficult to perform.
% FORM requires a limit state function to take the TSE of
% for me, my limit state function would be:
% g = E1 - 0.95
% pf = P(g<0)
% however, g must be a function of two or more randomly distrubuted normal
% variables. 
% My only variable is theta for each layer
% I have an equation explicitly solving for E1 for one layer (it's a big f)
% I could take the derivative of this function and then complete the
% hasofer lind algorithm to find the shortest distance.
% however, with only one degree of freedom, my function and derivative will
% be very complicated. 
% I would be wasting my time trying to develop functions for E1 for
% laminates of multiple layers
% is there another way to model all of this data in a less computationally
% expensive manner...?

%% Gaussian Process



