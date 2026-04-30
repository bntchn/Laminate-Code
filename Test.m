%% Recreating Results from Previous Publications
% [1] https://www.sciencedirect.com/science/article/pii/S1359836813005337
% [2] https://www.sciencedirect.com/science/article/pii/S1110016822008420#b0170

figure(1)
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

%% Monte Carlo Sampling of 1-5 layers
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

%% Visualizing Monte Carlo Data

pf = zeros(5,3);
failure_condition = 0.95;

for k = 1:3
    figure(k+1)
    for t = 1:layers
        E_normalized{k}(t,:) = sort(E_normalized{k}(t,:));
        loglog(E_normalized{k}(t,:)) 
        hold on

        data = E_normalized{k}(t,:);
        pf(t,k) = sum(data < failure_condition)/numel(data)*100;
    end
    hold off
end
xlabel('Monte-Carlo Samples'); ylabel('Normalized E_1')

%% Histograms of Monte Carlo Data

for k = 1:3 
    figure(k+4)
    for t = 1:layers
        subplot(5,1,t)
        histogram(E_normalized{k}(t,:),0:0.02:1)
        xlim([min(E_normalized{k}(t,:)),1])
    end
end


