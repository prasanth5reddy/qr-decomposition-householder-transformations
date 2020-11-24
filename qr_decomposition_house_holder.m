% set random seed for reproducibility
rng(1);
% initialize matrix
A = randi([1 30], 6, 5);
% uncomment below line to check with written example
% A = [-6 18 1; -3 9 3; 3 -29 -9];
disp('Given matrix A');
disp(A);

% compute QR decomposition using householder function
[Q, R] = qr_house_holder(A);
disp('QR decomposition using householder transformations')
disp('Q');
disp(Q);
disp('R');
disp(R);
% verify if QR properties are satisfied
if norm(Q * R - A) < 1e-10
    disp('Q * R = A is satifisied');
end
if norm(Q * Q' - eye(length(A))) < 1e-10
    disp(['Q * Q'' = I is satifisied', newline]);
end

% compute QR decomposition using built-in matlab function
% matlab function does not bother to make main diagonal entries of R +ve.
% Hence matrix elements signs might differ
[Q, R] = qr(A);
disp('QR decomposition using matlab function')
disp('Q');
disp(Q);
disp('R');
disp(R);
