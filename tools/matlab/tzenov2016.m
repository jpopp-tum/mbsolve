% tzenov2016  Displays result of tzenov2016 setup.
%             See https://doi.org/10.1364/OE.24.023232 for reference data.

% mbsolve: An open-source solver tool for the Maxwell-Bloch equations.
%
% Copyright (c) 2016, Computational Photonics Group, Technical University
% of Munich.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

clear all;
close all;

%  choose hdf5 file
[filename, folder] = uigetfile('../../*.hdf', 'Select result data');
f = fullfile(folder, filename);

% read global attributes
d_x = h5readatt(f, '/', 'gridpoint_size');
d_t = h5readatt(f, '/', 'timestep_size');
t_e = h5readatt(f, '/', 'sim_endtime');
L_x = h5readatt(f, '/', 'dev_length');

% complete grid
x = 0:d_x:L_x;
t = 0:d_t:t_e;

% read simulation data
e = h5read(f, '/e0/real');

% cut off transients (set start to 0 to use full trace)
start = 15e-9;
interval = t > start;
e = e(interval);
t = t(interval);

% plot field
figure;
plot(t, e(1, :));
xlabel('Time/s');
ylabel('Field at facet/a.u.');

% plot populations
d11 = h5read(f, '/d11/real');
d22 = h5read(f, '/d22/real');
d33 = h5read(f, '/d33/real');
d44 = h5read(f, '/d44/real');
d55 = h5read(f, '/d55/real');
d11 = d11(interval);
d22 = d22(interval);
d33 = d33(interval);
d44 = d44(interval);
d55 = d55(interval);
figure;
plot(t, d11, 'DisplayName', 'd11');
hold on;
plot(t, d22, 'DisplayName', 'd22');
plot(t, d33, 'DisplayName', 'd33');
plot(t, d44, 'DisplayName', 'd44');
plot(t, d55, 'DisplayName', 'd55');
legend show;

% plot intensity
I = abs(e).^2;
 
figure;
plot(t, I);
%xlim([1.8e-9, 1.84e-9]);
xlabel('Time/s');
ylabel('Field at facet/a.u.');

[Y, freqs] = calc_spectrum(e, d_t);

% spectrum
papersize = [ 15 6 ];
fig = figure('units', 'centimeters');
pos = get(gcf, 'pos');
semilogy(freqs/1e12, Y.^2/max(Y.^2), 'Color', [0, 101, 189]/255);
xlim([3.1, 4.1]);
xticks(3.1:0.1:4.1);
ylim([1e-4, 10]);
yticks(10.^[-4:1]);
xlabel('Frequency/THz');
ylabel('Power spectrum at facet/a.u.');

set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize', papersize);
print(fig, 'tzenov2016.pdf', '-dpdf', '-fillpage');

Y_I = calc_spectrum(I, d_t);

figure;
plot(freqs, 10*log10(Y_I));
xlim([4e9, 12e9]); 
xlabel('Frequency/Hz');
ylabel('Beatnote/a.u.');

return;

% apply filter
f_sep = 3.833e12;
filter_hi = (abs(f) > f_sep);
filter_lo = (abs(f) <= f_sep);

Y_hi = ifftshift(Y) .* filter_hi;
Y_lo = ifftshift(Y) .* filter_lo;

e_hi = real(fft(fftshift(Y_hi))) ./ win;
e_lo = real(fft(fftshift(Y_lo))) ./ win;

% plot parts
T_RT = 1/8.3e9;
rt_end = floor(t_e/T_RT);
rt_begin = rt_end - 3;

figure;
subplot(2, 1, 1);
plot(t/T_RT, e_lo .^2, '-', 'DisplayName', 'Low freq', 'Color', tum_colors(1, :))
xlabel('Time/T_{RT}');
ylabel('Intensity at facet/a.u.');
legend show;
xlim([rt_begin, rt_end]);
%hold on;
subplot(2, 1, 2);
plot(t/T_RT, e_hi .^2, '-', 'DisplayName', 'High freq', 'Color', tum_colors(2, :))
xlabel('Time/T_{RT}');
ylabel('Intensity at facet/a.u.');
legend show;
xlim([rt_begin, rt_end]);
