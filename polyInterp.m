f1 = @(x) sqrt(x.^3);
f2 = @(x) x .* sin(1 - x.^2);
f3 = @(x) 36 ./ (2 * x + 1).^3;
f4 = @(x) x .* exp(-x.^2);

limits = [
    0, 1;
    -1, 1;
    0, 1;
    1, 2
];

acc_vals = [
    integral(f1, 0, 1);
    integral(f2, -1, 1);
    integral(f3, 0, 1);
    integral(f4, 1, 2)
];


function result = midpoint_rule(f, a, b)
    result = (b - a) * f((a + b) / 2);
end

function result = trapezoid_rule(f, a, b, n)
    h = (b - a) / n;
    result = f(a) + f(b);
    for i = 1:n-1
        result = result + 2 * f(a + i * h);
    end
    result = result * h / 2;
end

function result = simpsons_rule(f, a, b, n)
    h = (b - a) / n;
    result = f(a) + f(b);
    for i = 1:2:n-1
        result = result + 4 * f(a + i * h);
    end
    for i = 2:2:n-2
        result = result + 2 * f(a + i * h);
    end
    result = result * h / 3;
end

function result = gaussian_quadrature(f, a, b)
    x1 = -1/sqrt(3);
    x2 = 1/sqrt(3);
    w1 = 1;
    w2 = 1;
    result = (b - a) / 2 * (w1 * f((b - a) / 2 * x1 + (a + b) / 2) + ...
                             w2 * f((b - a) / 2 * x2 + (a + b) / 2));
end

methods = {@midpoint_rule, @trapezoid_rule, @simpsons_rule, @gaussian_quadrature};
names = {'Midpoint', 'Composite Trapezoid', 'Simpson', 'Gaussian Quadrature'};

for i = 1:length(methods)
    fprintf('\n%s Rule Results:\n', names{i});
    for j = 1:4
        a = limits(j, 1);
        b = limits(j, 2);
        f = eval(['f', num2str(j)]);
        if i == 2 || i == 3
            result = methods{i}(f, a, b, 10);
        else
            result = methods{i}(f, a, b);
        end
        error = abs(result - acc_vals(j));
        fprintf('Integral (%g): Result = %.8f, Accurate = %.8f, Error = %.8f\n', ...
            j, result, acc_vals(j), error);
    end
end






f = @(x) 4 ./ (1 + x.^2);

a = 0;
b = 1;

function result = trapezoid_rule2(f, a, b, n)
    h = (b - a) / n;
    result = f(a) + f(b);
    for i = 1:n-1
        result = result + 2 * f(a + i * h);
    end
    result = result * h / 2;
end

function result = simpsons_rule2(f, a, b, n)
    h = (b - a) / n;
    result = f(a) + f(b);
    for i = 1:2:n-1
        result = result + 4 * f(a + i * h);
    end
    for i = 2:2:n-2
        result = result + 2 * f(a + i * h);
    end
    result = result * h / 3;
end

true_value = pi;

n_values = [4, 10];
fprintf('\n(a) Composite Trapezoidal Rule Results:\n');
for n = n_values
    approx_pi = trapezoid_rule2(f, a, b, n);
    error = abs(approx_pi - true_value);
    fprintf('n = %d: Approximate π = %.8f, Error = %.8f\n', n, approx_pi, error);
end

fprintf('\n(b) Simpson Rule Results:\n');
for n = n_values
    approx_pi = simpsons_rule2(f, a, b, n);
    error = abs(approx_pi - true_value);
    fprintf('n = %d: Approximate π = %.8f, Error = %.8f\n', n, approx_pi, error);
end
