%OFDM Example
%Example of how to use the OFDM class to generate a signal.
%
% Author: Chance Tarver
% Website: http://www.chancetarver.com
% July 2018;
%
% Modified:
%   Nov 2019. Add RFWebLab option.

%% ------------- BEGIN CODE --------------

use_rf_weblab = 1;  % if 0, we'll just add a little noise to the signal.

% Setup OFDM
params.nSubcarriers = 1200;
params.subcarrier_spacing = 15e3; % 15kHz subcarrier spacing
params.constellation = 'QPSK';
params.cp_length = 144; % Number of samples in cyclic prefix.
params.nSymbols = 1;
modulator = OFDM(params);

% Modulate and transmit 
[tx_data, tx_constellations] = modulator.use;
if use_rf_weblab
    dbm_power = -24;
    pa = webRF(dbm_power);
    % Upsample and send through RFWeblab and downsample
    [pa_out, pa_out_upsampled, pa_input_upsampled] = ...
        pa.transmit(tx_data, modulator.sampling_rate);
else
    pa_out = tx_data + 0.005 * randn(size(tx_data));
end
 
% Demod
rx_constellations = modulator.demod(pa_out);

% Error
error = rx_constellations - tx_constellations;
if norm(error)/length(error) < 0.01
    disp('Success! RX Data = TX Data');
else
    disp('Something may not be correct...');
end

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