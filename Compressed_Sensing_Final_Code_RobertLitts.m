%% Robert Litts, ECE 612 Digitial Signal Processing
%% Term Project, Spring 2020, Old Dominion University
%% Distributed Compressed Sensing of Correlated Sparse Signals

close all
clear
clc

%% SETUP: INITIALIZE SIGNAL CONSTANTS
rng(0);                 % set RNG seed
N = 256;                % length of signal
T =8;                  % number of non-zero peaks
M =46;                 % number of measurements to take (N < L)
%C = 3;                  %Number of users in network
SMNR=2000;                %Signal to measurement noise Ratio
plots = N/2;            %Number of measurements to plot

for C=1:9  
%Generate cells for each of "C" users to store data
x=cell(1,C) %Original signal
A=cell(1,C) %Measurement matrix
y=cell(1,C) %Transmitted signal (M measurements of original T-Sparse signal "x")
x0=cell(1,C) %initial guess at x using x=A'y
xp=cell(1,C) %output, reconstructed signal
error_1 = cell(1, C) %Error during reconstruction of original signal after CS
error = zeros(1, plots);

     
    
    

count = 1;
for M=1:plots %Vary the measurement values
support_set_voting = zeros(1,N) %Voting matrix to store democratic votes
common_support_set=zeros(1,T) %Common (Joint) Support set based on merging of individual support sets


%% GENERATE GAUSSIAN SPARSE SIGNALS, EACH W/ "T" NON-ZERO VALUES
for j=1:C %C users in network

    x{j} = zeros(N, 1); % original signal (T-sparse)
    peaks = randperm(N); % Generate signal with T randomly spread values
    peaks = peaks(1:T); %Determine the peak locations 
    support_set(j, 1:T) = peaks(1:T) %Partial Support set for each signal in network
    x{j}(peaks(1:T)) = randn(1,T) %Gaussian T-Sparse values set at peak locations
    amp(j) = 1.5*max(abs(x{j}));%Used to help plot signals
    x_none{j}=x{j};
    %x{j}=awgn(x{j}, SMNR); %Add noise to signal, 20dB SNR, this will be the received signal
    
    % Generate Sensing matrix MxN (Measurements = M, Signal size = N,
    % additional measurements for simulated additional users in network
    A{j} = sqrt(1/(M+(C-1)))*randn(M+(C-1),N);
    
    %Taking M random measurements of original signal 
    y{j} = A{j}*x{j};   
end


%% CONSENSUS VOTING
for p = 1:C %C Users in network
    [~, t0] = sort(abs(x{p}), 'descend')
    for q=1:T
    support_set_voting(1,t0(q)) = support_set_voting(1,t0(q)) + 1 %check every value within each signal's individual support set and create democratic voting
    end
end


%% FUSION
adder = 1
for r = 1:N  %Create Common Support set from individual support sets based on democratic voting
    if support_set_voting(1,r) >= 2
        common_support_set(1,adder) = r
         adder = adder + 1
    end
end

%% EXPAND COMMON SUPPORT SET 
if length(find(common_support_set)) < T %If less than T non-zero values, use T values
    t1 = common_support_set(1:T)
else 
    [~, t0] = sort(abs(common_support_set), 'descend') %if more than T values, lexographically sort 
    %t1 = t0(1:T)
    t1 = common_support_set(1:T)
end
common_support_set(1:T) = t1;


%% MERGE COMMON SUPPORT SET WITH INDIVIDUAL SIGNALS

for j=1:C
x0{j} = A{j}.'*y{j}; %Initial Guess for signal

check=0;
if sum(common_support_set)>0  
for k=1:N   
for m=1:T
    if common_support_set(1,m) ==k && common_support_set(1,m) >0
        x0{j}(common_support_set(1,m)) = 100*x0{j}(common_support_set(1,m));
        check = 1;
    elseif common_support_set(1,m) == 0
        check =1;
    end      
end

    if check ==0
        x0{j}(k) = x0{j}(k)*.01;
    end
end
end
 
xp{j} = l1eq_pd(x0{j}, A{j}, [], y{j}); %[2]
error_1{j} =  mean((xp{j}-x{j}).^2);

if j==C
error(count) = sqrt(mean((x{j}-xp{j}).^2));
count = count + 1;
end
    
end

end




u = (1:plots)
figure(4)
if C==1
plot(u/N, smoothdata(error(1:plots), 'gaussian',50), 'r--', 'LineWidth', 2, 'DisplayName', ['C', num2str(C)']); title('Error of Recovered Signal Compared to Original Signal'); ylim([0 .3])
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal$$', 'Interpreter','latex','FontSize',10);
xlabel('$$ \alpha(M/N)$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
hold on

elseif C==2
plot(u/N, smoothdata(error(1:plots), 'gaussian',50), 'b--', 'DisplayName', ['C', num2str(C)']); title('Error of Recovered Signal Compared to Original Signal'); ylim([0 .3])
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal$$', 'Interpreter','latex','FontSize',10);
xlabel('$$ \alpha(M/N)$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
hold on

elseif C==3
plot(u/N, smoothdata(error(1:plots), 'gaussian',50), 'm-.', 'DisplayName', ['C', num2str(C)']); title('Error of Recovered Signal Compared to Original Signal'); ylim([0 .3])
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal$$', 'Interpreter','latex','FontSize',10);
xlabel('$$ \alpha(M/N)$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
grid on
hold on

elseif C==4
plot(u/N, smoothdata(error(1:plots), 'gaussian',50), ':g.', 'DisplayName', ['C', num2str(C)']); title('Error of Recovered Signal Compared to Original Signal'); ylim([0 .3])
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal$$', 'Interpreter','latex','FontSize',10);
xlabel('$$ \alpha(M/N)$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
grid on
hold on


elseif C==9
plot(u/N, smoothdata(error(1:plots), 'gaussian',50), 'k-','DisplayName', ['C', num2str(C)']); title('Error of Recovered Signal Compared to Original Signal'); ylim([0 .3])
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal$$', 'Interpreter','latex','FontSize',10);
xlabel('$$ \alpha(M/N)$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
grid on
hold on

end
end
    
figure(1) %Plot original signal
subplot(2,1,1); plot(x_none{3}); title(['$$ Original\ Gaussian\ Sparse\ Signal\ x_1:\ T=$$',num2str(T), ' Sparse'], 'Interpreter','latex','FontSize',10); xlim([1 N]); ylim([-amp(1) amp(1)]);
xlabel('$$ N$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Amplitude$$', 'Interpreter','latex','FontSize',10);
subplot(2,1,2); plot(x{3}); title(['$$ Original\ Gaussian\ Sparse\ Signal\ With\ Noise\ x_1:\ T=$$',num2str(T), ' Sparse,\ SNR=20dB\ '], 'Interpreter','latex','FontSize',10); xlim([1 N]); ylim([-amp(1) amp(1)]);
xlabel('$$N$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Amplitude$$', 'Interpreter','latex','FontSize',10);
grid on

figure(2) %Plot random measurements of original signal baased on sensing matrix A
plot(y{1}); title(['$$ "M"\ Random\ Measurements\ of\ Gaussian\ Sparse\ Signal\ x_1:\ T=$$',num2str(T), ' Sparse'], 'Interpreter','latex','FontSize',10); ylim([-amp(1) amp(1)]);
xlabel('$$M\ (Measurement)$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Amplitude$$', 'Interpreter','latex','FontSize',10);
grid on

figure(3); plot(xp{3}, 'DisplayName', 'Recovered Signal'); xlim([1 N]); ylim([-amp(1) amp(1)]);
hold on
plot(x{3}, '--', 'DisplayName', 'Original Signal')
legend show
title('$$Recovered\ Signal$$', 'Interpreter','latex','FontSize',10);
xlabel('$$ N$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Amplitude$$', 'Interpreter','latex','FontSize',10);
grid on
       

%% TEST BELOW FOR VARYING SMNR W/ FIXED ALPHA=.18
%% SETUP: INITIALIZE SIGNAL CONSTANTS
rng(0);                 % set RNG seed
N = 256;                % length of signal
T =8;                  % number of non-zero peaks
M =46;                 % constant number of measurements, alpha =.18
plots = 30;            %Number of measurements to plot

for C=1:9  
%Generate cells for each of "C" users to store data
x=cell(1,C) %Original signal
A=cell(1,C) %Measurement matrix
y=cell(1,C) %Transmitted signal (M measurements of original T-Sparse signal "x")
x0=cell(1,C) %initial guess at x using x=A'y
xp=cell(1,C) %output, reconstructed signal
error_1 = cell(1, C) %Error during reconstruction of original signal after CS
error = zeros(1, plots);

     
    
    

count = 1;
for SMNR=1:30 %Vary the measurement values
support_set_voting = zeros(1,N) %Voting matrix to store democratic votes
common_support_set=zeros(1,T) %Common (Joint) Support set based on merging of individual support sets


%% GENERATE GAUSSIAN SPARSE SIGNALS, EACH W/ "T" NON-ZERO VALUES
for j=1:C %C users in network

    x{j} = zeros(N, 1); % original signal (T-sparse)
    peaks = randperm(N); % Generate signal with T randomly spread values
    peaks = peaks(1:T); %Determine the peak locations 
    support_set(j, 1:T) = peaks(1:T) %Partial Support set for each signal in network
    x{j}(peaks(1:T)) = randn(1,T) %Gaussian T-Sparse values set at peak locations
    amp(j) = 1.5*max(abs(x{j}));%Used to help plot signals
    x_none{j}=x{j};
    x{j}=awgn(x{j}, SMNR); %Add noise to signal, varying the SNR
    
    % Generate Sensing matrix MxN (Measurements = M, Signal size = N,
    % additional measurements for simulated additional users in network
    A{j} = sqrt(1/(M+(C-1)))*randn(M+(C-1),N);
    
    %Taking M random measurements of original signal 
    y{j} = A{j}*x{j};   
end


%% CONSENSUS VOTING
for p = 1:C %C Users in network
    [~, t0] = sort(abs(x{p}), 'descend')
    for q=1:T
    support_set_voting(1,t0(q)) = support_set_voting(1,t0(q)) + 1 %check every value within each signal's individual support set and create democratic voting
    end
end


%% FUSION
adder = 1
for r = 1:N  %Create Common Support set from individual support sets based on democratic voting
    if support_set_voting(1,r) >= 2
        common_support_set(1,adder) = r
         adder = adder + 1
    end
end

%% EXPAND COMMON SUPPORT SET 
if length(find(common_support_set)) < T %If less than T non-zero values, use T values
    t1 = common_support_set(1:T)
else 
    [~, t0] = sort(abs(common_support_set), 'descend') %if more than T values, lexographically sort 
    %t1 = t0(1:T)
    t1 = common_support_set(1:T)
end
common_support_set(1:T) = t1;


%% MERGE COMMON SUPPORT SET WITH INDIVIDUAL SIGNALS

for j=1:C
x0{j} = A{j}.'*y{j}; %Initial Guess for signal

check=0;
if sum(common_support_set)>0  
for k=1:N   
for m=1:T
    if common_support_set(1,m) ==k && common_support_set(1,m) >0
        x0{j}(common_support_set(1,m)) = 100*x0{j}(common_support_set(1,m));
        check = 1;
    elseif common_support_set(1,m) == 0
        check =1;
    end      
end

    if check ==0
        x0{j}(k) = x0{j}(k)*.01;
    end
end
end
 
xp{j} = l1eq_pd(x0{j}, A{j}, [], y{j}); %[2]
error_1{j} =  mean((xp{j}-x{j}).^2);

%combine = cat(3,error_1{:});
%total_error(count) = mean(combine,3);
%error(count) = sqrt(mean((x{j}-xp{j}).^2));
%error1(count) = sqrt(mean((x{j}).^2));
%error(r) = sqrt(mean((xp{j}-x{1:j}).^2));
%SRER(count) = error1(count)/error(count)
if j==C
error(count) = sqrt(mean((x{j}-xp{j}).^2));
count = count + 1;
end
    
end

end




u = (1:plots)
figure(5)
if C==1
plot(u, smoothdata(error(1:plots), 'gaussian',50), 'r--', 'LineWidth', 2, 'DisplayName', ['C', num2str(C)']); title('Error of Recovered Signal Compared to Original Signal'); 
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal, \alpha=.18$$', 'Interpreter','latex','FontSize',10);
xlabel('$$SMNR(dB)$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
hold on

elseif C==2
plot(u, smoothdata(error(1:plots), 'gaussian',50), 'b--', 'DisplayName', ['C', num2str(C)']); title('Error of Recovered Signal Compared to Original Signal');
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal, \alpha=.18$$', 'Interpreter','latex','FontSize',10);
xlabel('$$SMNR(dB)$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
hold on

elseif C==3
plot(u, smoothdata(error(1:plots), 'gaussian',50), '-.*', 'DisplayName', ['C', num2str(C)']); title('Error of Recovered Signal Compared to Original Signal'); 
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal, \alpha=.18$$', 'Interpreter','latex','FontSize',10);
xlabel('$$ SMNR(dB) $$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
grid on
hold on

elseif C==4
plot(u, smoothdata(error(1:plots), 'gaussian',50), ':b.', 'DisplayName', ['C', num2str(C)']); title('Error of Recovered Signal Compared to Original Signal'); 
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal, \alpha=.18$$', 'Interpreter','latex','FontSize',10);
xlabel('$$SMNR(dB) $$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
grid on
hold on


elseif C==9
plot(u, smoothdata(error(1:plots), 'gaussian',50), '-gd', 'DisplayName', ['C', num2str(C)']); title('Error of Recovered Signal Compared to Original Signal'); 
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal, \alpha=.18$$', 'Interpreter','latex','FontSize',10);
xlabel('$$SMNR(dB) $$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
grid on
hold on
end

end


%% TEST BELOW FOR WATTS-STROGATZ NETWORK MODEL
%% SETUP: INITIALIZE SIGNAL CONSTANTS
rng(0);                 % set RNG seed
N = 256;                % length of signal
T =8;                  % number of non-zero peaks
M =46;                 % constant number of measurements, alpha =.18
plots = N/8;            %Number of measurements to plot
SMNR=20;                %Set SMNR at 20dB

%% Watts Strogatz Large Random Network
q = 3;
p=.3 %probability of being rewired to uniformly chosen random node
C = 100; %100 nodes in ranom network
Watts_Matrix = rand(1,C);
Watts_nodes = find(Watts_Matrix>(1-p)); %Get nodes which are >(1-p) probability

for C=1:length(Watts_nodes)  
%Generate cells for each of "C" users to store data
x=cell(1,C) %Original signal
A=cell(1,C) %Measurement matrix
y=cell(1,C) %Transmitted signal (M measurements of original T-Sparse signal "x")
x0=cell(1,C) %initial guess at x using x=A'y
xp=cell(1,C) %output, reconstructed signal
error_1 = cell(1, C) %Error during reconstruction of original signal after CS
error = zeros(1, plots);


        

     
    
    

count = 1;
for M=1:plots %Vary the measurement values
support_set_voting = zeros(1,N) %Voting matrix to store democratic votes
common_support_set=zeros(1,T) %Common (Joint) Support set based on merging of individual support sets


%% GENERATE GAUSSIAN SPARSE SIGNALS, EACH W/ "T" NON-ZERO VALUES
for j=1:C %C users in network

    x{j} = zeros(N, 1); % original signal (T-sparse)
    peaks = randperm(N); % Generate signal with T randomly spread values
    peaks = peaks(1:T); %Determine the peak locations 
    support_set(j, 1:T) = peaks(1:T) %Partial Support set for each signal in network
    x{j}(peaks(1:T)) = randn(1,T) %Gaussian T-Sparse values set at peak locations
    amp(j) = 1.5*max(abs(x{j}));%Used to help plot signals
    x_none{j}=x{j};
    x{j}=awgn(x{j}, SMNR); %Add noise to signal, varying the SNR
    
    % Generate Sensing matrix MxN (Measurements = M, Signal size = N,
    % additional measurements for simulated additional users in network
    A{j} = sqrt(1/(M+(C-1)))*randn(M+(C-1),N);
    
    %Taking M random measurements of original signal 
    y{j} = A{j}*x{j};   
end


%% CONSENSUS VOTING
for p = 1:C %C Users in network
    [~, t0] = sort(abs(x{p}), 'descend')
    for q=1:T
    support_set_voting(1,t0(q)) = support_set_voting(1,t0(q)) + 1 %check every value within each signal's individual support set and create democratic voting
    end
end


%% FUSION
adder = 1
for r = 1:N  %Create Common Support set from individual support sets based on democratic voting
    if support_set_voting(1,r) >= 2
        common_support_set(1,adder) = r
         adder = adder + 1
    end
end

%% EXPAND COMMON SUPPORT SET 
if length(find(common_support_set)) < T %If less than T non-zero values, use T values
    t1 = common_support_set(1:T)
else 
    [~, t0] = sort(abs(common_support_set), 'descend') %if more than T values, lexographically sort 
    %t1 = t0(1:T)
    t1 = common_support_set(1:T)
end
common_support_set(1:T) = t1;


%% MERGE COMMON SUPPORT SET WITH INDIVIDUAL SIGNALS

for j=1:C
x0{j} = A{j}.'*y{j}; %Initial Guess for signal

check=0;
if sum(common_support_set)>0  
for k=1:N   
for m=1:T
    if common_support_set(1,m) ==k && common_support_set(1,m) >0
        x0{j}(common_support_set(1,m)) = 100*x0{j}(common_support_set(1,m));
        check = 1;
    elseif common_support_set(1,m) == 0
        check =1;
    end      
end

    if check ==0
        x0{j}(k) = x0{j}(k)*.01;
    end
end
end
 
xp{j} = l1eq_pd(x0{j}, A{j}, [], y{j}); %[2]
error_1{j} =  mean((xp{j}-x{j}).^2);

%combine = cat(3,error_1{:});
%total_error(count) = mean(combine,3);
%error(count) = sqrt(mean((x{j}-xp{j}).^2));
%error1(count) = sqrt(mean((x{j}).^2));
%error(r) = sqrt(mean((xp{j}-x{1:j}).^2));
%SRER(count) = error1(count)/error(count)
if j==C
error(count) = sqrt(mean((x{j}-xp{j}).^2));
count = count + 1;
end
    
end

end




u = (1:plots)
figure(5)
if C==1
plot(u/N, smoothdata(error(1:plots), 'gaussian',50), 'r--', 'LineWidth', 2, 'DisplayName', 'C1'); title('Error of Recovered Signal Compared to Original Signal'); 
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal$$', 'Interpreter','latex','FontSize',10);
xlabel('$$ \alpha(M/N)$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
hold on


elseif C==length(Watts_nodes)
plot(u/N, smoothdata(error(1:plots), 'gaussian',50), '-gd', 'DisplayName', ['Watts-Strogatz Model']); title('Error of Recovered Signal Compared to Original Signal'); 
title('$$Error\ of\ Recovered\ Gaussian\ Signal\ Compared\ to\ Original\ Signal$$', 'Interpreter','latex','FontSize',10);
xlabel('$$ \alpha(M/N)$$', 'Interpreter','latex','FontSize',10);
ylabel('$$Error$$', 'Interpreter','latex','FontSize',10);
legend show
grid on
hold on
end

end


%% REFERENCES
%% [1] D. Sundman, S. Chatterjee, and M. Skoglund, “Design and Analysis of a Greedy Pursuit for Distributed Compressed Sensing,” IEEE Transactions on Signal Processing, vol. 64, no. 11, pp. 2803–2818, 2016. 
%% [2] https://statweb.stanford.edu/~candes/software/l1magic/index.html#links    
   

