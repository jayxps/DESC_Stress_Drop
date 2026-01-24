function [doy, year] = julian2doy(jd)
%  JULIAN2DOY Calculates the day of the year (and Gregorian year) of a
%   given Julian date.
%   DOY = JULIAN2DOY( JD ) converts one or more dates into the day of the
%   year. Input JD can be an M-by-1 or 1-by-M vector containing M julian
%   dates. JULIAN2DOY returns a vector (same dimensions) of M day-of-the-year
%   numbers.
%
%   Examples:
%
%   Calculate doy of the year of 2454088.5 (December 19, 2006), giving 353:
%	   [doy, year] = julian2doy(2454088.5)
%
%   Calculate day of the year of row vector of two dates, giving [352, 353]:
%	   doy = julian2doy([2454088.4 2454089.5])
%
%   Calculate day of the year of column vector of same dates, [353; 353]:
%	   doy = julian2doy([2454088.9; 2454089.0])
%
%   Notes:
%    - Conversion modified from julian2greg.m by Gabriel Ruiz Mtz,
%      Jun-2006.
%    - Function also returns year, because alone doy is ambiguous.
%
%   See also MJULIAN2DOY, JULIANDATE, MJULIANDATE.
%   Copyright 2016 Jacco Geul <jacco@geul.net>, redistribution and use,
%      with or without modification, are permitted exclusively under the
%      terms of the Modified BSD license. See LICENSE.txt file.
% Find out the year
I = floor( jd + 0.5);
A = floor( ( I- 1867216.25 ) / 36524.25 );
a4 = floor( A / 4 );
B = I;
B(I>=2299160) = I + 1 + A - a4;
C = B + 1524;
D = floor( ( C - 122.1 ) / 365.25 );
E = floor( 365.25 * D );
G = floor( ( C - E ) / 30.6001 );
month = G - 13;
month(G <= 13.5) = G - 1;
year = D - 4715;
year(month > 2.5) = D - 4716;
% Create two zero and one matrices of same size.
m0 = zeros(size(year));
m1 = ones(size(year));
% Compute the day of the year
doy = floor(jd - juliandate(year,m1,m1,m0,m0,m0)) + 1;

