% Finds value of k0 by calculating the harmonic partial sum H(k) where the
% sum ceases to increase any further due to limitations of single precision
% The harmonic series theoretically diverges but due to the limitations of
% single precision and rounding errors of floating-point. The sum
% H(k) = sum(1/i for i=1 to k) which stops increasing when 1/k becomes
% smaller than epsilon. This machine epsilon lies around 1.192 * 10^-7 and
% the terms 1/k become too small to be represented one k exceeds a certain
% number and the harmonic sum stops increasing.
% We can use the approximation H(k) ~ log(k) for large k and knowing 1k ~
% machine epsilon, we estimate that the harmonic sum stops increasing
% around k0 ~ 8.39 * 10^6





prev = single(0);
curr = single(0);

k = 1;

while true
    curr = curr + single(1/k);

    if curr == prev
        k0 = k - 1;  
        break;
    end
    prev = curr;
    k = k + 1;
end

fprintf('The sum stops increasing at k0 = %d\n', k0);