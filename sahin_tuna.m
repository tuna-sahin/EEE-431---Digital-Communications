%% %Q1-b-1
clf;
u = rand(1,100)+rand(1,100);
histogram(u-1,65)
title('100 Experiments')
%% %Q1-b-2
clf;
u = rand(1,100000) + rand(1,100000);
histogram(u-1,65)
title('100000 Experiments')
x=u-1;
%% Q1-d
clf;
Ns = [8 16 64];
for i = 1:3
    subplot(1,3,i)
    [boundaries,midpoints] = get_uniform_quantizer(-1,1,Ns(i));
    x_hat = quantize(boundaries,midpoints,x);
    [mse_error,sqnr] = get_sqnr(x,x_hat);
    histogram(x_hat,129)
    title(strcat("MSE Error: ",num2str(mse_error,'%0.2f'),"     SQNR: ",num2str(sqnr,'%0.2f')))
end
%% Q1-e 
% Calculate the area of each trapezoid
clf;
Ns = [8 16 64];
for j = 1:3
    subplot(1,3,j)
    [boundaries,midpoints] = get_uniform_quantizer(-1,1,Ns(j));
    probabilities = zeros(size(midpoints));

    for i = 1:length(boundaries)-1
        probabilities(i) = abs(boundaries(i+1)-boundaries(i))*(1-abs(boundaries(i))+1-abs(boundaries(i+1)))/2;
    end
    stem(midpoints,probabilities,'filled','MarkerSize',3)
    title(strcat('PMF of Q(X) with N=',num2str(N)))
end
%% Q1-f
[boundaries_8,midpoints_8] = lloyd_max(-1,1,8);
[boundaries_16,midpoints_16] = lloyd_max(-1,1,16);
[boundaries_64,midpoints_64] = lloyd_max(-1,1,64);
%% Obtain the plot for Quantization Functions
clf;
x_2 = linspace(-1,1,1000);
subplot(1,3,1)
fineplot(x_2,quantize(boundaries_8,midpoints_8,x_2),'Q(X)','x','Q(x)',[-1 1],[-1 1],'off',[400 400],'','-');
subplot(1,3,2)
fineplot(x_2,quantize(boundaries_16,midpoints_16,x_2),'Q(X)','x','Q(x)',[-1 1],[-1 1],'off',[400 400],'','-');
subplot(1,3,3)
fineplot(x_2,quantize(boundaries_64,midpoints_64,x_2),'Q(X)','x','Q(x)',[-1 1],[-1 1],'off',[400 400],'','-');
%% Q1-g Obtain histograms and sqnrs 
clf;
for i = 1:3
    subplot(1,3,i)
    switch i
        case 1
            boundaries = boundaries_8;
            midpoints = midpoints_8;
        case 2
            boundaries = boundaries_16;
            midpoints = midpoints_16;
        case 3
            boundaries = boundaries_64;
            midpoints = midpoints_64;
        otherwise
            disp("error")
    end
    x_hat = quantize(boundaries,midpoints,x);
    [mse_error,sqnr] = get_sqnr(x,x_hat);
    histogram(x_hat,129)
    title(strcat("MSE Error: ",num2str(mse_error,'%0.2f'),"     SQNR: ",num2str(sqnr,'%0.2f')))
end
%% Q-1h
clf;
Ns = [8 16 64];
for j = 1:3
    switch j
        case 1
            boundaries = boundaries_8;
            midpoints = midpoints_8;
        case 2
            boundaries = boundaries_16;
            midpoints = midpoints_16;
        case 3
            boundaries = boundaries_64;
            midpoints = midpoints_64;
        otherwise
            disp("error")
    end
    subplot(1,3,j)
    probabilities = zeros(size(midpoints));

    for i = 1:length(boundaries)-1
        probabilities(i) = abs(boundaries(i+1)-boundaries(i))*(1-abs(boundaries(i))+1-abs(boundaries(i+1)))/2;
    end
    stem(midpoints,probabilities,'filled','MarkerSize',3)
    title(strcat('PMF of Q(X) with N=',num2str(N)))
end
%% Q2-a
%% Read First Image
I = imread("foliage.tif");
imshow(I)
%% Read Second Image
I = imread("testimage.tif");
imshow(I,[])
%% Plot Image Histogram
clf;
imhist(I)
%% Obtain and Plot Polynomial Approximation
clf;
histogram_values = imhist(I);
x1 = linspace(0,255,length(histogram_values));
p = polyfit(x1,histogram_values,4);
hold on
approx = polyval(p,x1);
plotted = (approx .* (approx > 0)) + 10^-12;
fineplot(x1,plotted/sum(plotted),'PDF Approximation of Fourth-Order','Pixel Intensity (I)','F_I(i)',[-25 255],[-0.001 0.01],'off',[400 100],'off','-')
grid on

%% Q2-b
arrayified_I = reshape(I,[1 numel(I)]); %sorry for bad variable names lol
quantization_steps = [2 4 8 16 32];
for i = 1:length(quantization_steps)
    [boundaries,midpoints] = get_uniform_quantizer(0,255,quantization_steps(i));
    quantized = quantize(boundaries,midpoints,arrayified_I);
    subplot(2,3,i)
    I2 = reshape(quantized,size(I));
    title(strcat('Quantization intervals = ',num2str(quantization_steps(i))))
    imshow(I2,[]);
end
%% Calculate sqnrs and mse
arrayified_I = reshape(I,[1 numel(I)]); %sorry for bad variable names lol
quantization_steps = [2 4 8 16 32];
mse = [];
sqnr = [];
for i = 1:length(quantization_steps)
    [boundaries,midpoints] = get_uniform_quantizer(0,255,quantization_steps(i));
    quantized = quantize(boundaries,midpoints,arrayified_I);
    mse = [mse sum((double(arrayified_I)-quantized).^2)];
    sqnr = [sqnr 10*log10(sum(double(arrayified_I).^2)/sum((double(arrayified_I)-quantized).^2))];
end
sqnr 
mse
clf;
finestem(1:5,sqnr,'SQNR as a Function of Quantization Steps','N','SQNR',[-1 6],[-5 45],'off',[400 400],'','-');
xticks([1 2 3 4 5])
xticklabels({'2','4','8','16','32'})


%% Q2-c
quantization_steps = [2 4 8 16 32];
for i = 1:length(quantization_steps)
    [boundaries,midpoints] = get_uniform_quantizer(0,255,quantization_steps(i));
    quantized = quantize(boundaries,midpoints,arrayified_I);
    subplot(2,3,i)
    I2 = reshape(quantized,size(I));
    title(strcat('Quantization intervals = ',num2str(quantization_steps(i))))
    imshow(trunc(fftshift(20*log10(10^-12+abs(fft2(I2))))),[])
end
    %% Q2-d
clf;
[boundaries,midpoints] = get_uniform_quantizer(0,255,8);
I_ft = fftshift(fft2(reshape(quantize(boundaries,midpoints,arrayified_I),size(I))));
a = 0.1;
sizes_ft = size(I_ft);
mask_x = abs(linspace(-1,1,sizes_ft(1))) < a;
mask_y = abs(linspace(-1,1,sizes_ft(2))) < a;
mask = (mask_x.' * mask_y);
subplot(1,2,1)
imshow(ifft2(fftshift(I_ft.*mask)),[])
[boundaries,midpoints] = get_uniform_quantizer(0,255,8);
I_ft = fftshift(fft2(reshape(arrayified_I,size(I))));
a = 0.1;
sizes_ft = size(I_ft);
mask_x = abs(linspace(-1,1,sizes_ft(1))) < a;
mask_y = abs(linspace(-1,1,sizes_ft(2))) < a;
mask = (mask_x.' * mask_y);
subplot(1,2,2)
imshow(ifft2(fftshift(I_ft.*mask)),[])

%% Functions
function y = trunc(x) % rounds everything below -40dB up to -40dB. Maps -40,0 to 0,1
    y = ((x+40).*(x+40>0)./40);
end
function [boundaries,midpoints] = lloyd_max(lower_bound,upper_bound,N)
    boundaries = linspace(lower_bound,upper_bound,N+1);
    boundaries = boundaries(2:end-1); %eliminate -1 1 so they dont messup the algorithm
    midpoints = zeros(1,length(boundaries)+1);
    
    size(boundaries)
    for iter = 1:10000
        boundaries = [-1 boundaries 1];
        for i = 1:length(boundaries)-1
            
            samples = linspace(boundaries(i),boundaries(i+1),1000);
            pdf = 1-abs(samples); % change this depending on the pdf
            pdf = pdf/sum(pdf); %normalise for the conditional pdf
            
            midpoints(i) = sum(samples.*pdf); %take the expectation
        end
        boundaries = boundaries(2:end-1);
        
        for i = 1:length(midpoints)-1
            boundaries(i) = (midpoints(i)+midpoints(i+1))/2;
        end
    end
    boundaries = [-1 boundaries 1];
end
function [mse, sqnr] = get_sqnr(x,x_hat)
    power = sum(x.^2);
    mse = sum((x-x_hat).^2);
    sqnr = 10*log10(power/mse);
end
function [boundaries,midpoints] = get_uniform_quantizer(lower_bound,upper_bound,N)
    boundaries = linspace(lower_bound,upper_bound,N+1);
    midpoints = (boundaries(1:end-1) + boundaries (2:end))/2;
end
function x_hat = quantize(boundaries,midpoints,x)

x_hat = zeros(size(x));
for i = 1:length(x)
    for j = 1:length(boundaries)-1
        if x(i) >= boundaries(j) && x(i) <= boundaries(j+1)
            x_hat(i) = midpoints(j);
        end
    end
end
end