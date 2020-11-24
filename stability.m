% set random seed for reproducibility
rng(1);

% initialize variables for plotting
x = 1:100;
norm_cgs = zeros(size(x));
norm_mgs = zeros(size(x));
norm_hh = zeros(size(x));

for n=x
    % set up a ill-conditioned matrix using auto generated SVD
    singularVec = 0.5 * 0.5.^(0:n-1);
    sigma = diag(singularVec);
    U = randn(n, n);
    V = randn(n, n);
    A = U * sigma * V';

    % compute QR using Classical Gram-Schmidt
    [Q, R] = cgs(A);
    norm_cgs(n) = norm(Q * Q' - eye(n));

    % compute QR using Modified Gram-Schmidt
    [Q, R] = mgs(A);
    norm_mgs(n) = norm(Q * Q' - eye(n));

    % compute QR using House Holder reflections
    [Q, R] = qr_house_holder(A);
    norm_hh(n) = norm(Q * Q' - eye(n));
end

hold on;
plot(x, norm_cgs, '-r', 'LineWidth', 1.5);
plot(x, norm_mgs, '-b', 'LineWidth', 1.5);
plot(x, norm_hh, 'color', [0 0.5 0], 'LineWidth', 1.5);
legend('Classical Gram-Schimdt', 'Modified Gram-Schimdt', 'HouseHolder');
xlabel('Size of Matrix');
ylabel('Orthoganility Error');
title('HouseHolder Stability');



