%OFDM Example
%Example of how to use the OFDM class to generate a signal.
%
% Author: Chance Tarver
% Website: http://www.chancetarver.com
% July 2018;
%
% Modified:
%   Jan 2020. Add windowing bw symbols, absolute cyclic prefix, and evm
%   Nov 2019. Add RFWebLab option.

%% ------------- BEGIN CODE --------------

use_rf_weblab = 1;  % if 0, we'll just add a little noise to the signal.

% Setup OFDM
params.n_subcarriers = 1200;
params.subcarrier_spacing = 15e3; % 15kHz subcarrier spacing
params.constellation = 'QPSK';
params.n_symbols = 1;
params.use_windowing = true;
modulator = OFDM(params);

% Modulate and transmit 
[tx_data, ~] = modulator.use;
if use_rf_weblab
    dbm_power = -30;
    pa = webRF(dbm_power);
    % Upsample and send through RFWeblab and downsample
    [pa_out, pa_out_upsampled, pa_input_upsampled] = ...
        pa.transmit(tx_data, modulator.sampling_rate);
else
    pa_out = tx_data + 0.005 * randn(size(tx_data));
end
 
% Demod
rx_constellations = modulator.demod(pa_out);

% Error Vector Magnitude
evm = modulator.calculate_evm(rx_constellations);
fprintf(' EVM = %d %%\n', evm);

modulator.plot_rx_errors(rx_constellations);

figure
plot(real(pa_out))
grid on
title('Time Domain View')
xlabel('Sample')
ylabel('Magnitude')

figure
plot(rx_constellations, 'o')
hold on
plot(tx_constellations, 'o')
title('RX and TX Constellations')