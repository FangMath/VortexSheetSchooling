function pot = freeByBodyHil(zetaf, zetabody, nu, Para)
    N = length(zetabody) - 1;
    M = length(zetaf);

    %val(1) = nu(1)./sqrt(1-Para.s(2).^2); % get values
    %val(N+1) = nu(N+1)./sqrt(1-Para.s(N).^2); % get values
    %val(2:N) = nu(2:N)./sqrt(1-Para.s(2:N).^2); % get values

    cbody = zetabody(N/2+1); % mid body point
    shat = (zetabody(end) - zetabody(1))/2;

    target = (zetaf - cbody)/shat;
if M >1
    target(M) = target(M-1);
    else
    target(M) = target(M) + 1e-9;
end

    pot = hilbert_sig_cheb(nu, target);
    %pot = hilbert_sig(val, target, 1);
    pot = pot/shat;

    pot = pot/(2*pi*1i); % get velocity
    pot = conj(pot);
end
