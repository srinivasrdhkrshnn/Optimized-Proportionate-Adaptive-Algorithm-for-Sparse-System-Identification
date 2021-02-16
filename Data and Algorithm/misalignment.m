
%% EE6110 Project: Author - K.R.SRINIVAS EE18B136

% Code A : The below program allows you plot and compare the misalignment
% for different algorithms.

clear all
close all

load IR_sparse.mat;  % loads the echo path model
load speech.mat ;    % loads the far end signal; can use if synthetic signal not reqd

%% Genaration of synthetic speech sequence
 
a = 0.8 ;
Nr = [sqrt(1-a^2)] ;
Dr = [1 -a] ;
N = 600000 ;                 % sequence length
u = randn(N,1);             
far_end = filter(Nr,Dr,u);   % input sequence of unit variance (AR(1) process)
var_n = 0.01 ;               % Near end Background noise
N = length(far_end);

%% Impulse Response of the system

h1 = IR_sparse ;
h2 = circshift(h1,50) ;      % shifted impulse response to model change in echo path model

%% Echo/Desired signal

echo_1 = filter(h1,1,far_end(1:N/2)) ;
echo_2 = filter(h2,1,far_end((N/2)+1:end)) ;
echo = [echo_1;echo_2] + sqrt(var_n)*randn(N,1) ;  % desired signal (echo + near end)

%% Initializations
M = 512 ;                    % Filter-Tap Length

% IPNLMS Specifications mu = 0.1
mu = 0.1 ;
alpha = 0 ;
delta = 9.8314e-04 ;
epsilon = 0.01 ;

w0 = zeros(M,1) ;            % Weight vector
u0 = zeros(M,1) ;            % regressor vector
m0 = zeros(N,1) ;            % misalignment vector


disp('Please wait for a while...')

%% IPNLMS Algorithm

for i=1:N
    u0 = [far_end(i);u0(1:M-1)];   % Regressor vector update
    e0(i) = echo(i) - u0'*w0;      % apriori error
  
    for s = 1:M
        k(s) = (1-alpha)/2*M + (1+alpha)*norm(w0(s),1)/(2*norm(w0,1)+delta) ; % proportionate step-size implementation      
    end
    
    Q = diag(k) ;                  % Step-Size update matrix
    
    w0 = w0 + (mu*e0(i)*Q*u0)/(u0'*Q*u0 + epsilon) ;    % weight update rule

    if i <= N/2                                         % misalignment
         m0(i) = 20*log10(norm(h1-w0)/norm(h1)) ;
    else
         m0(i) = 20*log10(norm(h2-w0)/norm(h2)) ;
    end 
    
    if mod(i,5000)==0
         i
    end
end

% IPNLMS Specifications  mu = 1
mu = 1 ;
alpha = 0 ;
delta = 9.8314e-04 ;
epsilon = 0.01 ;

w3 = zeros(M,1) ;            % Weight vector
u3 = zeros(M,1) ;            % regressor vector
m3 = zeros(N,1) ;            % misalignment vector


disp('Please wait for a while...')
%% IPNLMS Algorithm

for i=1:N
    u3 = [far_end(i);u3(1:M-1)];   % Regressor vector update
    e3(i) = echo(i) - u3'*w3;      % apriori error
  
    for s = 1:M
        k(s) = (1-alpha)/2*M + (1+alpha)*norm(w3(s),1)/(2*norm(w3,1)+delta) ; % proportionate step-size implementation      
    end
    
    Q = diag(k) ;                  % Step-Size update matrix
    
    w3 = w3 + (mu*e3(i)*Q*u3)/(u3'*Q*u3 + epsilon) ;   % weight update rule

    if i <= N/2                                        % misalignment
         m3(i) = 20*log10(norm(h1-w3)/norm(h1)) ;
    else
         m3(i) = 20*log10(norm(h2-w3)/norm(h2)) ;
    end 
    
    if mod(i,5000)==0
         i
    end
end

% OPLMS Specifications

w2 = zeros(M,1) ;            % Weight vector
u2 = zeros(M,1) ;            % regressor vector
m2 = zeros(N,1) ;            % misalignment matrix
m = 1e-2 ;
var_w = 0 ;                  % process noise
gamma = ones(M,1);
I_l = ones(M,1)  ;

%% OPLMS Algorithm

for i = 1:N
   u2 = [far_end(i) ; u2(1:M-1)] ;         % Regressor vector update
   e2(i) = echo(i)-u2'*w2 ;                % apriori error
   
   var_x = (u2'*u2)/M  ;                   % variance of input signal
   q = M /(m + M*var_w);
   mu_new = 1/(q*var_n + var_x*M) ;
   w2_old = w2 ;
   w2 = w2 + q*mu_new*(gamma.*u2)*e2(i) ;  % weight vector update rule
   p = gamma ;
   gamma = gamma + var_w*I_l + var_x*(1-2*(q^2)*(mu_new^2))*(gamma.*gamma) ;
   r = max(q*gamma) ;
   gamma = (1/r)*gamma ;
   m = m + M*var_w - q*mu_new*var_x*norm(p)^2 ;
   
   var_w = (1/M)*norm(w2-w2_old)^2 ;      % process noise variance update
  
   
   if i <= N/2                            % misalignment
         m2(i) = 20*log10(norm(h1-w2)/norm(h1)) ;
   else
         m2(i) = 20*log10(norm(h2-w2)/norm(h2)) ;
   end 
    
   if mod(i,5000)==0
        i
   end
end

plot(1:N,m0); hold on
plot(1:N,m3); hold on
plot(1:N,m2); 
xlabel('sample sequence');
ylabel('Misalignment (dB)');
legend('IPNLMS : mu = 0.1','IPNLMS : mu = 1','OPLMS');
grid

% figure
% subplot(311)
% plot(1:N,far_end);
% title('Far-end signal');
% axis([0 N -10 10]);
% 
% 
% subplot(312)
% plot(1:N,echo);
% title('Echo signal');
% axis([0 N -10 10]);
% 
% 
% subplot(313)
% plot(1:N,e3);
% title('Error signal');
% axis([0 N -10 10]);

