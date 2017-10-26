%Probability for poisson distribution 
%lambda (lam) should be less than 10 
function val = poisson(n,lam)
    global poissonBackup;
    key = n*10 +lam;
    key_as_char = int2str(key);
    if ~poissonBackup.isKey(key_as_char)
        poissonBackup(key_as_char) = exp(-lam)*(lam^n)/factorial(n);
    end
    val = poissonBackup(key_as_char);
end