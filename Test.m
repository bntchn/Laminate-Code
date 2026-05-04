%% Recreating Results from Previous Publications
% [1] https://www.sciencedirect.com/science/article/pii/S1359836813005337
% [2] https://www.sciencedirect.com/science/article/pii/S1110016822008420#b0170

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
xlabel('Orientation, \theta (deg)'); ylabel('Normalized E_1')

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

%% Monte Carlo Data Parameters

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

%% Visualizing Monte Carlo Data
% this is not a cdf, maybe clarify that with him
% how could this be turned into a cdf?

for k = 1:3
    figure
    for t = 1:layers
        E_normalized{k}(t,:) = sort(E_normalized{k}(t,:));
        plot(E_normalized{k}(t,:)) 
        hold on
    end
    hold off
end
xlabel('Monte-Carlo Samples'); ylabel('Normalized E_1')

%% Histograms of Monte Carlo Data
% Talk to Mirzendehdel about this

% would it be worth it to create graphs and characterize my orientaion
% data?

% what should my bin count be? I have the formula, but should the bin edges
% start and end where my data is? or should it start and end at 0 and 1.
% My "distribution" can look differently based on what my bin edges are.
% what is the best way to map my bin edges?

% what should my y-axis be? should it be pdf probability divided by bin
% width (graphs with values over 100)? or should it be the probability 
% that a sample falls in that bin (graphs with 0 to 1)?

% I can use histogram function to make probability, pdf, and cdf graphs

bin_count = 5*10*round((1 + 3.3*log10(sample))/10);
edges = linspace(0,1,bin_count+1);

for k = 1:3 
    figure
    for t = 1:layers
        subplot(layers,1,t)
        histogram(E_normalized{k}(t,:),edges,'Normalization','probability')
    end
    ylabel('Probability (%)'); xlabel('Normalized E_1')
end

%% Matching to a Distribution
% plotting the most extreme distribution, data = E_normalized{3}(1,:)
% these plots are matching pdf, not probability

omega = 1;          % upper bound
lambda = 0.3;       % scale parameter (controls spread)
k = 0.5;            % shape parameter
x = 0:0.001:0.999;

z = (omega - x) / lambda;

f = (k / lambda) .* z.^(k-1) .* exp(-(z.^k));

figure
plot(x, f, 'LineWidth', 2)
xlabel('x')
ylabel('f(x)')
title('Reverse Weibull Distribution (bounded at \omega)')
grid on
hold on
histogram(E_normalized{3}(1,:),edges,'Normalization','pdf')

% chi-squared test
% compare pdf

% KS test
% compare cdf

%% FORM

% I do not want to derive the explicit relationship of E1 for the purposes
% of doing a FORM analysis.
% P( g < 0 ) , g = E1(normalized, theta1, theta2, theta3...) - 0.95
% that sounds terrible. I already have an explicit relationship for E1 of
% one layer as a function of theta, given by Ref [2].
% I would propose to Dr. Mirzendehdel that I can do this instead. I would
% compare the time it takes to do a proper monte carlo simulation of 1
% layer vs HL algorithm of one layer. 
% This would be for the purpose of showing I know what FORM is.


