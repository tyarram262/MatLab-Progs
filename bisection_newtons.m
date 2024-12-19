equations = {
    @(x) x^3 - 2*x - 5, @(x) 3*x^2 - 2, [2, 4], 4;
    @(x) exp(-x) - x, @(x) -exp(-x) - 1, [0, 2], 2;
    @(x) x * sin(x), @(x) sin(x) + x*cos(x), [-1, 1], 1;
    @(x) x^3 - 3*x^2 + 3*x - 1, @(x) 3*x^2 - 6*x + 3, [-1, 1], 0;
};

tolerance = 1e-6;
max_iter = 50;

figure;
titles = {'x^3 - 2x - 5 = 0', 'e^{-x} - x = 0', 'x sin(x) = 0', 'x^3 - 3x^2 + 3x - 1 = 0'};

for i = 1:4
    f = equations{i, 1};
    df = equations{i, 2};
    interval = equations{i, 3};
    x0 = equations{i, 4};
    
    [bisection_errors, ~] = bisection_method(f, interval(1), interval(2), tolerance, max_iter);
    
    [newton_errors, ~] = newtons_method(f, df, x0, tolerance, max_iter);
    
    subplot(2, 2, i);
    hold on;
    plot(1:length(bisection_errors), bisection_errors, '-o', 'LineWidth', 1.5, 'DisplayName', 'Bisection');
    plot(1:length(newton_errors), newton_errors, '-x', 'LineWidth', 1.5, 'DisplayName', 'Newton');
    title(titles{i});
    xlabel('Iterations');
    ylabel('Error');
    legend('Location', 'northeast');
    grid on;
end


function [errors, root] = bisection_method(f, a, b, tol, max_iter)
    errors = [];
    for k = 1:max_iter
        c = (a + b) / 2;
        errors(k) = abs(f(c));
        if abs(f(c)) < tol
            root = c;
            return;
        elseif f(a) * f(c) < 0
            b = c;
        else
            a = c;
        end
    end
    root = c;
end

f = @(x) x^5 - x^3 - 4*x;
df = @(x) 5*x^4 - 3*x^2 - 4;

x0_a = 1;
x0_b = 2;
tolerance = 1e-6;
max_iter = 50;

[errors_a, root_a, iter_a] = newtons_method(f, df, x0_a, tolerance, max_iter);

[errors_b, root_b, iter_b] = newtons_method(f, df, x0_b, tolerance, max_iter);

fprintf('Part (a): Initial value x0 = 1\n');
fprintf('Root: %.6f, Iterations: %d\n', root_a, iter_a);
fprintf('Part (b): Initial value x0 = 2\n');
fprintf('Root: %.6f, Iterations: %d\n', root_b, iter_b);

figure;
hold on;
plot(1:length(errors_a), errors_a, '-o', 'LineWidth', 1.5, 'DisplayName', 'x0 = 1');
plot(1:length(errors_b), errors_b, '-x', 'LineWidth', 1.5, 'DisplayName', 'x0 = 2');
title('Error Convergence of Newton’s Method');
xlabel('Iterations');
ylabel('Error (|f(x)|)');
legend('Location', 'northeast');
grid on;

function [errors, root, num_iter] = newtons_method(f, df, x0, tol, max_iter)
    errors = [];
    x = x0;
    for k = 1:max_iter
        fx = f(x);
        errors(k) = abs(fx);
        if abs(fx) < tol
            root = x;
            num_iter = k;
            return;
        end
        dfx = df(x);
        if dfx == 0
            error('Derivative became zero.');
        end
        x = x - fx / dfx;
    end
    root = x;
    num_iter = max_iter;
end


% When x0 = 1 the Newton's method exhibits oscillatory behavior due to the
% derivative being to small or alternating in sign close to x0 = 1 which
% may cause instability
%This could either cause it to fail to converge or take a significant
%amount of iterations.

%Some unusual observations would be that since x0=1 lives in a region where
%f'(x) approaches zero this causes large corrections in x and move the
%guess further from the root. This may cause hte iterations to diverge or
%stagnate depending on f'(x).

%Starting with x0 = 2 Newtons method will converge to a root and the
%derivative is well-defined and nonzero in the region in the region around
%x0=2 ensuring rapid convergence.

%The root near x=2 is found successfully after a few iterations though.