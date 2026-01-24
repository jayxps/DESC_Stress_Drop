function [fc, fc2, o0, Es, delsig, asig, r, G, err, dspecfit, iflag, ftest, etest, fcb1, fcb2, falloff] ...
    = sourcepara_vs(dspec, f, nfrq, imethod, moment, vs, iwave)

global falloffin
% Now adding double-corner frequency fitting

% Medium parameters (using variable shear velocity vs)
% beta = 4030;
% beta = 3464;
beta = vs;
pho = 2700;
mu = pho * beta^2;
kfact = 0.26;

% Get set of source parameters based on wave type
if (iwave == 1)
    % P-wave parameters
    % kfact = 0.38;
    % fact = 16./7 * (0.38*beta).^3;  % this is for P-wave, k = 0.32..S k = 0.21~ 0.2766
    fact = (0.42 * beta)^3;  % to be consistent with fortran code
    alpha = 1.732 * beta;
    scale = 2. / (15 * pi * pho * alpha^5) * 14.4;
else
    % S-wave parameters
    kfact = 0.26;
    % fact = 16./7 * (0.27*beta).^3;
    fact = (0.2766 * beta)^3.;  % this is for P-wave, k = 0.32..S k = 0.21~ 0.2766
    % fact = (0.42*beta)^3.;
    scale = 1.07 / (5 * pi * pho * beta^5);
end

% Optimization options
options = optimset('TolX', 1e-10, 'MaxFunEvals', 1500, 'MaxIter', 1500);

% If this is not resampled, the subroutine will resample this to equal spacing in log frequency domain
fmin = min(f(nfrq));
fmax = max(f(nfrq));
dfx = diff(f);
dfxmean = abs(dfx - mean(dfx));

if (max(dfxmean) < 1e-3)  % linear spaced
    t = logspace(log10(fmin), log10(fmax), 100);
    y = spline(f(nfrq), dspec(nfrq), t);
else  % already log spaced
    t = f(nfrq);
    y = dspec(nfrq);
end

% Method setup for fitting dspec
im = imethod;
npar = 2;
falloff = falloffin;

if (imethod == 3)
    im = 5;
    npar = 3;
    falloff = falloffin;  % boatwright model
end
if (imethod == 2)
    im = 1;
    npar = 3;
    falloff = falloffin;  % brune model
end
if (imethod == 4)
    im = 6;
    npar = 4;
end
if (imethod == 8)
    im = 7;
    npar = 4;
end
if (imethod == 9)
    im = 10;
    npar = 4;
end
if (imethod == 10 || imethod == 11)
    npar = 3;
end

% Set initial parameters based on npar
nfx = find(f >= 2 & f <= 4);
x0 = mean(dspec(nfx));

if (npar == 3)
    pstart = [x0 10 50];
elseif (npar == 2)
    pstart = [x0 10];
elseif (npar == 4)
    pstart = [x0 10 50 2];
end

if (imethod == 6 | imethod == 7)
    npar = 3;
    pstart = [x0 10 2];
end

if (imethod == 10)
    pstart = [x0 10 2];  % initial fall off is 2
end

if (imethod == 11)
    pstart = [x0 2 20];
end

% Optimization: standard fminsearch or grid search for double corner
if (imethod ~= 11)
    [p_method3, ~, iflag] = fminsearch(@(x)method(x, t, y, imethod), pstart, options);
elseif (imethod == 11)
    % Grid search for the double corner method
    errmin = 99999;
    fc = 0.1:0.1:100;
    for ifc1 = 1:length(fc)
        for ifc2 = ifc1+1:length(fc)
            fc1 = fc(ifc1);
            fc2 = fc(ifc2);
            upred = -1/2 * log10(1 + (t/fc1).^2) - 1/2 * log10(1 + (t/fc2).^2);
            udiff = y - upred;
            omega0 = mean(udiff);
            ndiff = udiff - omega0;
            err = norm(ndiff);
            if (err < errmin)
                errmin = err;
                p_method3(2) = fc1;
                p_method3(3) = fc2;
                p_method3(1) = omega0;
            end
        end
    end
    % u = p(1) - 1/2*log10(1+(t/p(2)).^2) - 1/2* log10(1+(t/p(3)).^2);
end

% Calculate error and fit
err = method(p_method3, t, y, imethod);
u3 = fitmethod(p_method3, f, imethod);

% Compute variance reduction (commented out)
% nf = length(f);
% nfx = nf*2-1;
% dx(1:nf) = u3;
% for i = nf+1:nfx
%     dx(i) = u3(nfx-i+1);
% end

% Extract fitted parameters
dspecfit = u3;
o0 = p_method3(1);
fc = abs(p_method3(2));  % corner frequency for target event

if (npar == 3)
    fc2 = abs(p_method3(3));  % corner frequency for egf
elseif (npar == 2)
    fc2 = 1000;
elseif (npar == 4)
    fc2 = abs(p_method3(3));
    falloff = p_method3(4);
end

if (imethod == 6 | imethod == 7)
    fc2 = 1000;
end

if (imethod == 10)
    fc2 = 1000;
    falloff = p_method3(3);
end

% Test corner frequency uncertainty
fc1 = fc - 5;
if (fc1 <= 0)
    fc1 = 0.05;
end
ftest = linspace(fc1, fc+5, 50);

for i = 1:length(ftest)
    % if (npar == 3)
    %     p_test = [o0 fc2];
    % else
    %     p_test = [o0 ];
    % end
    p_test = p_method3;

    [p_t, fval, iflag] = fminsearch(@(x)testmethod(x, ftest(i), t, y, imethod), p_test, options);

    etest(i) = fval;
end

% Find corner frequency bounds
% minetest = min(etest);
% if (minetest ~= 0) etest = etest/minetest; end
nlow = find(etest <= 1.05*min(etest));
% plot(ftest, etest)

if (~isempty(nlow))
    fcb1 = ftest(nlow(1));
    fcb2 = ftest(nlow(length(nlow)));
else
    fcb1 = -1;
    fcb2 = -1;
end

% Fill in parameter array
p(1) = log10(moment);
p(2) = p_method3(2);

% Fill in fall off rate
if (imethod == 4)
    p(3) = p_method3(4);
    falloff = p(3);
end
if (imethod == 8)
    p(3) = p_method3(4);
    falloff = p(3);
end
if (imethod == 6)
    falloff = p_method3(3);
    p(3) = falloff;
end
if (imethod == 7)
    falloff = p_method3(3);
    p(3) = falloff;
end
if (imethod == 9)
    falloff = p_method3(4);
    p(3) = falloff;
end
if (imethod == 10)
    falloff = p_method3(3);
    p(3) = falloff;
end
if (imethod == 11)
    p(3) = p_method3(3);
end

% Calculate energy for each station using velocity spectra
dspecx = dspec - o0 + p(1);

% Method is to integrate between fmin and fmax using spectra,
% then extrapolate to before and after
% nfrqx = find(f >= 2 & f <= 25);
nfrqx = nfrq;

[Es, fex, dspecex, dspecfitx] = energy(f, dspecx, nfrqx, p, im);

% Calculate stress drop and related parameters
delsig = fc^3 * moment / fact / 1e6;  % in MPa
Es = Es .* scale;
asig = mu * Es(1) / moment / 1e6;  % asig in MPa
r = kfact * beta / fc;
A = pi * r^2;
D = moment / (mu * A);
G = 0.5 * (delsig - 2*asig) .* D * 1e6;


% Nested function: test method with fixed corner frequency
function err = testmethod(p, p2, t, y, im)
global falloffin

% Model descriptions
% brune family: 1, 2, 7, 8, 11
% boatwright family: 3, 4, 5, 6, 9, 10

if (im == 1)
    u = p(1) - log10(1 + (t/p2).^falloffin);
elseif (im == 2)
    u = p(1) - log10(1 + (t/p2).^falloffin) + log10(1 + (t/p(3)).^falloffin);
elseif (im == 3)
    u = p(1) - log10(1 + (t/p2).^(2*falloffin))/2 + log10(1 + (t/p(3)).^(2*falloffin))/2;
elseif (im == 4)
    u = p(1) - 1./p(4)*log10(1 + (t/p2).^(2*p(4))) + 1./p(4)*log10(1 + (t/p(3)).^(2*p(4)));
elseif (im == 5)
    u = p(1) - log10((1 + (t/p2).^(2*falloffin)).^(1/2));
elseif (im == 6)
    u = p(1) - log10((1 + (t/p2).^(2*p(3))).^(1/p(3)));
elseif (im == 7)
    u = p(1) - log10(1 + (t/p2).^p(3));
elseif (im == 8)
    u = p(1) - log10(1 + (t/p2).^p(4)) + log10(1 + (t/p(3)).^p(4));  % brune ratio with n float
elseif (im == 9)
    u = p(1) - 1./2*log10(1 + (t/p2).^(2*p(4))) + 1./2*log10(1 + (t/p(3)).^(2*p(4))); % boatwright ratio with n float
elseif (im == 10)
    u = p(1) - 1./2*log10(1 + (t/p2).^(2*p(3)));  % boatwright with n float
elseif (im == 11)
    u = p(1) - 1/2*log10(1 + (t/p(2)).^2) - 1/2*log10(1 + (t/p(3)).^2);
    % u = p(1) - log10(1 + (t/p(2)) + (t/p(3)).^2);
    % the double-corner frequency model p(2) would be the original Brune's corner, p(3) would be the new corner
    % Ideally, for a perfect hauskull model, the w-2 falloff corner correspond to the rise time
    % the w-1 falloff corner correspond to the rupture duration
    % for larger events, w-1 corner is much smaller due to long rupture duration
end

err = norm(u - y);


% Nested function: main method for fitting
function err = method(p, t, y, im)
global falloffin

% Model descriptions
% brune family: 1, 2, 7, 8, 11
% boatwright family: 3, 4, 5, 6, 9, 10

if (im == 1)
    u = p(1) - log10(1 + (t/p(2)).^(falloffin));  % pure brune model
elseif (im == 2)
    u = p(1) - log10(1 + (t/p(2)).^(falloffin)) + log10(1 + (t/p(3)).^(falloffin));  % brune model with ratio
elseif (im == 3)
    u = p(1) - log10(1 + (t/p(2)).^(2*falloffin))/2 + log10(1 + (t/p(3)).^(2*falloffin))/2;  % boatwright with ratio
elseif (im == 4)
    u = p(1) - 1./p(4)*log10(1 + (t/p(2)).^(2*p(4))) + 1./p(4)*log10(1 + (t/p(3)).^(2*p(4)));  % boatwright ratio with gamma float
elseif (im == 5)
    u = p(1) - log10((1 + (t/p(2)).^(2*falloffin)).^(1/2));  % boatwright (no ratio)
elseif (im == 6)
    u = p(1) - log10((1 + (t/p(2)).^(2*p(3))).^(1/p(3)));  % boatwright with gamma float (no ratio)
elseif (im == 7)
    u = p(1) - log10(1 + (t/p(2)).^p(3));  % brune with n float (no ratio)
elseif (im == 8)
    u = p(1) - log10(1 + (t/p(2)).^p(4)) + log10(1 + (t/p(3)).^p(4));  % brune ratio with n float
elseif (im == 9)
    u = p(1) - 1./2*log10(1 + (t/p(2)).^(2*p(4))) + 1./2*log10(1 + (t/p(3)).^(2*p(4)));  % boatwright ratio with n float
elseif (im == 10)
    u = p(1) - 1./2*log10(1 + (t/p(2)).^(2*p(3)));  % boatwright with n float
elseif (im == 11)
    u = p(1) - 1/2*log10(1 + (t/p(2)).^2) - 1/2*log10(1 + (t/p(3)).^2);  % the double-corner frequency model (p(2) is second, p(3) is first)
end

err = norm(u - y);


% Nested function: fit method for prediction
function u = fitmethod(p, t, im)
global falloffin

if (im == 1)
    u = p(1) - log10(1 + (t/p(2)).^(falloffin));
elseif (im == 2)
    u = p(1) - log10(1 + (t/p(2)).^(falloffin)) + log10(1 + (t/p(3)).^(falloffin));
elseif (im == 3)
    u = p(1) - log10(1 + (t/p(2)).^(2*falloffin))/2 + log10(1 + (t/p(3)).^(2*falloffin))/2;
elseif (im == 4)
    u = p(1) - 1./p(4)*log10(1 + (t/p(2)).^(2*p(4))) + 1./p(4)*log10(1 + (t/p(3)).^(2*p(4)));
elseif (im == 5)
    u = p(1) - log10((1 + (t/p(2)).^(2*falloffin)).^(1/2));
elseif (im == 6)
    u = p(1) - log10((1 + (t/p(2)).^(2*p(3))).^(1/p(3)));
elseif (im == 7)
    u = p(1) - log10(1 + (t/p(2)).^p(3));
elseif (im == 8)
    u = p(1) - log10(1 + (t/p(2)).^p(4)) + log10(1 + (t/p(3)).^p(4));  % brune ratio with n float
elseif (im == 9)
    u = p(1) - 1./2*log10(1 + (t/p(2)).^(2*p(4))) + 1./2*log10(1 + (t/p(3)).^(2*p(4)));  % boatwright ratio with n float
elseif (im == 10)
    u = p(1) - 1./2*log10(1 + (t/p(2)).^(2*p(3)));  % boatwright ratio with n float
elseif (im == 11)
    u = p(1) - 1/2*log10(1 + (t/p(2)).^2) - 1/2*log10(1 + (t/p(3)).^2);
end


% Nested function: estimate energy from spectra
% Method: for frequency in nfrq, direct integration
% then extrapolate spectra using method 3
function [Esx, fex, dspecex, dspecfit] = energy(f, dspec, nfrq, p, im)

% Direct integration over measured frequency range
w0 = f(nfrq) * 2 * pi;
m0 = (w0 .* 10.^dspec(nfrq)).^2;
Es0 = trapz(f(nfrq), m0);

% Extrapolate below measured range
df = 1e-3;
fmin = 1e-10;
fmax = 2000;

f1 = fmin:df:f(nfrq(1))-df;
u1 = fitmethod(p, f1, im);
w1 = f1 * 2 * pi;
m1 = (w1 .* 10.^u1).^2;

% Extrapolate above measured range
f2 = f(nfrq(length(nfrq)))+df:df:fmax;
u2 = fitmethod(p, f2, im);
w2 = f2 * 2 * pi;
m2 = (w2 .* 10.^u2).^2;

% Full modeled range
f3 = fmin:df:fmax;
u3 = fitmethod(p, f3, im);
w3 = f3 * 2 * pi;
m3 = (w3 .* 10.^u3).^2;

% Integrate extrapolated portions
Es1 = trapz(f1, m1);
Es2 = trapz(f2, m2);
Es3 = trapz(f3, m3);

% Combine frequency and spectra arrays
fex = [f1, f(nfrq), f2];
dspecex = [u1, dspec(nfrq), u2];
dspecfit = fitmethod(p, fex, im);

% Total energy calculation
Es = Es0 + Es1 + Es2;

% Esx: [total energy, modeled energy, measured energy, extrapolated low, extrapolated high]
Esx = [Es Es3 Es0 Es1 Es2];
