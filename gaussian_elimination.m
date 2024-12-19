function x = gaussian_elimination(A, b)
    [n, m] = size(A);
    if n ~= m
        error('Matrix A is square');
    end    
    if det(A) == 0
        error('Matrix A is singular');
    end
    
    pivot = 1:n;    
    for k = 1:n-1
        [~, max_index] = max(abs(A(k:n, k)));
        max_index = max_index + k - 1;
        
        if max_index ~= k
            A([k, max_index], :) = A([max_index, k], :);
            b([k, max_index]) = b([max_index, k]);
            pivot([k, max_index]) = pivot([max_index, k]);
        end        
        for i = k+1:n
            factor = A(i, k) / A(k, k);
            A(i, k:n) = A(i, k:n) - factor * A(k, k:n);
            b(i) = b(i) - factor * b(k);
        end
    end
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (b(i) - A(i, i+1:n) * x(i+1:n)) / A(i, i);
    end
end
