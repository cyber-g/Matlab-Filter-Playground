% filter_experiment.m
% This script designs a simple FIR low-pass filter. 
% It generates a complex Gaussian test signal, applies the filter, and compares 
% the output using two methods: 
% 1) A banded Toeplitz matrix of the signal.
% 2) A banded Toeplitz matrix of the filter coefficients.
% The script visualizes the results and computes the Mean Squared Error (MSE) 
% for both methods to evaluate their accuracy.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: filter
% Author: Germain Pham
% email: cygerpham@free.fr
% April 2025; Last revision: 
% 
% Special thanks to Denis Gilbert, who provided the original code for "M-file Header Template"
% https://fr.mathworks.com/matlabcentral/fileexchange/4908-m-file-header-template
% which is adapted here for the purpose of this example.
% 
%------------- BEGIN CODE --------------

% Design a simple FIR low pass filter
% Define the filter specifications
Fs = 1000; % Sampling frequency
Fc = 100;  % Cut-off frequency
N  = 50;   % Filter order

% Generate the low pass FIR filter coefficients using Hamming window
b = fir1(N, Fc/(Fs/2), hamming(N+1));
b = b(:); % Ensure b is a column vector

% Generate a test signal
Nsamples = 2000; % Number of samples
Nrandomize = 10; % Number of realizations
x = randn(Nsamples, Nrandomize) + 1j * randn(Nsamples, Nrandomize); % Complex Gaussian noise
x = x ./ (ones(Nsamples,1)*sqrt(var(x))); % Normalize the signal

% Apply the filter to the test signal
y_ref_filter = filter(b, 1, x); % Filter the signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Method 1: Using the banded Toeplitz SIGNAL matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the Tril-toeplitz matrix of the signal
Tx = zeros(Nsamples,length(b),Nrandomize);
for nrand = 1:Nrandomize
    Tx_tmp = conj(tril(toeplitz(x(:,nrand)))); % Toeplitz matrix
    Tx(:,:,nrand) = Tx_tmp(:,1:length(b));     % Keep only the first N columns corresponding to the filter length
end

% Compute the filter output using matrix multiplication
y_matrix_x = squeeze(pagemtimes(Tx,repmat(b(:),[1 1 Nrandomize]))); % Matrix multiplication

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Method 2: Using the banded Toeplitz FILTER matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tb = tril(toeplitz([b;zeros(length(x)-length(b),1)])); % Toeplitz matrix for filter coefficients

% Compute the filter output using matrix multiplication
y_matrix_b = Tb * x; % Matrix multiplication


% Plots
nexttile
plot(real([y_ref_filter(:,1) y_matrix_x(:,1) y_matrix_b(:,1)]))
nexttile
plot(imag([y_ref_filter(:,1) y_matrix_x(:,1) y_matrix_b(:,1)]))
legend('Reference Filter', 'Matrix Method (x)', 'Matrix Method (b)')
title('Filter Output Comparison')
xlabel('Sample Index')
ylabel('Amplitude')
grid on


% Errors 
MSE_x = sqrt(mean(abs(y_ref_filter - y_matrix_x).^2)); % Mean Squared Error for x
MSE_b = sqrt(mean(abs(y_ref_filter - y_matrix_b).^2)); % Mean Squared Error for b

disp(['MSE for banded Toeplitz SIGNAL: ', num2str(MSE_x)]);
disp(['MSE for banded Toeplitz FILTER: ', num2str(MSE_b)]);

%------------- END OF CODE --------------
