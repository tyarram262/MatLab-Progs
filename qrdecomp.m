function [Q,R] = qr(A)
    [m,n] = size(A);
    if m<n
        error('Matrix must have at least as many rows as cols');
    end    
    if m == n && det(A) == 0
        error('Matrix is not full rank');
    end
    Q = eye(m);
    R = A;
    for k = 1:n
        x = R(k:m, k);
        if norm(x)== 0
            error('Matrix not full rank');
        end
        e1 = zeros(length(x), 1);
        e1(1) = norm(x);
        vec = x - e1;
        vec = vec / norm(vec);
        R(k:m, k:n) = R(k:m, k:n) - 2 * (vec * (vec' * R(k:m, k:n)));

        Qnew = eye(m);
        Qnew(k:m, k:m) = Qnew(k:m, k:m) - 2 * (vec * vec');
        Q = Q * Qnew;
    end
    disp('Q =');
    disp(Q);
    disp('R =');
    disp(R);
end
t = [-1, -0.75, -0.5, 0, 0.25, 0.5, 0.75]';
y = [1.0000, 0.8125, 0.7500, 1.0000, 1.3125, 1.7500, 2.3125]';
A = [t.^2, t, ones(length(t), 1)];
[Q, R] = qr(A);
x = R \ (Q' * y);
x2 = A \ y;
residual = norm(y - A*x);
res2 = norm(y - A*x2);
disp('QR decomp Residual norm:');
disp(residual);
disp('Matlab Residual norm:');
disp(res2);
x = [2.2, 2.6, 3.4, 4.0]';
y = [65, 61, 54, 50]';
X = log(x);
Y = log(y);
A = [X, ones(length(X), 1)];
[Q, R] = qr(A);
coeffs_qr = R \ (Q' * Y);
b_qr = coeffs_qr(1);
A_qr = coeffs_qr(2);
a_qr = exp(A_qr);
disp(['a = ', num2str(a_qr)]);
disp(['b = ', num2str(b_qr)]);
% These are the 2 data points that would minimize the 2 norm of the
% residuals and we find that they are 92.4813 for a and -0.44137 for b