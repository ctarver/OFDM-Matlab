classdef OFDM < handle
    %OFDM Class for an OFDM signal.
    %
    %	Author:	Chance Tarver (2018)
    %		tarver.chance@gmail.com
    %
    properties
        n_subcarriers
        subcarrier_spacing
        constellation
        symbol_alphabet
        cp_length
        cp_legnth_start_of_slot
        n_symbols
        fft_size
        sampling_rate
        use_rrc_window
        window_length
        rrc_taps
        last_tx_time_d % previous transmission in time domain
        last_tx_freq_d % previous transmission in freq domain
    end
    
    properties (Constant, Hidden)
        constellation_library = {'QPSK','16QAM','64QAM'};
        constellaton_order=[4 16 64];
        n_active_scs = [72 180 300 600 900 1200];
        window_lengths=[4 6 4 6 8 8];  % https://www.mathworks.com/help/lte/ref/lteofdmmodulate.html#bugx3kl-1_head
        fft_sizes = [128 256 512 1024 2048]
        subcarrier_spacings = [15 30 60 120 240];
        cp_lengths_us_normal = [4.69 2.34 1.17 0.57 0.29]; % length of cp in microseconds for each numerology
    end
    
    methods
        function obj = OFDM(params)
            %OFDM Construct an instance of this class. Will create an OFDM
            %signal in the frequency and time domain. Will also upsample for PA
            %
            % Args:
            %    params: struct with all the settings.
            %
            %	Author:	Chance Tarver (2018)
            %		tarver.chance@gmail.com
            %
            
            obj.n_subcarriers =  params.n_subcarriers;
            obj.subcarrier_spacing = params.subcarrier_spacing;
            obj.constellation = params.constellation;
            obj.n_symbols = params.n_symbols;
            
            obj.fft_size = 2^ceil(log2(obj.n_subcarriers));
            obj.sampling_rate = obj.subcarrier_spacing * obj.fft_size;
            period = obj.sampling_rate^-1;
            cp_dictionary = containers.Map(obj.subcarrier_spacings, ...
                obj.cp_lengths_us_normal);
            cp_length_us = cp_dictionary(obj.subcarrier_spacing/1000)*1e-6;
            obj.cp_length = round(cp_length_us / period);
            fprintf(' Cyclic prefix = %d samples', obj.cp_length);
            obj.symbol_alphabet = obj.QAM_Alphabet(obj.constellation);
            obj.use_rrc_window = params.use_windowing;
            
            if params.use_windowing
                window_length_dictionary = containers.Map(obj.n_active_scs, ...
                    obj.window_lengths);
                try
                    obj.window_length = window_length_dictionary(obj.n_subcarriers);
                catch
                    obj.window_length = 8;
                end
                obj.generate_rrc();
            else
                obj.window_length = 0;
            end
        end
        
        
        function [out, fd_symbols] = use(obj)
            %use. Use the modulator.
            
            % Create random symbols on the constellation
            fd_symbols = zeros(obj.n_subcarriers, obj.n_symbols);
            fd_symbols = obj.symbol_alphabet(ceil(length(obj.symbol_alphabet) ...
                * rand(size(fd_symbols))));
            
            % iterate over symbols.
            td_symbols = zeros(obj.fft_size + obj.cp_length + obj.window_length, obj.n_symbols);
            for i = 1:obj.n_symbols
                td_waveform = obj.frequency_to_time_domain(fd_symbols(:, i));
                cp_td_waveform = obj.add_cyclic_prefic(td_waveform);
                td_symbols(:, i) = obj.add_windows(cp_td_waveform);
            end
            out = obj.create_full_waveform(td_symbols);
            obj.last_tx_freq_d = fd_symbols;
            obj.last_tx_time_d = out;            
        end
        
        function out = create_full_waveform(obj, in)
            N = obj.window_length;
            K = obj.n_symbols;
            samp_per_sym = obj.fft_size+obj.cp_length;
            out = zeros((samp_per_sym)*obj.n_symbols, 1);
            
            % first symbol is special
            out(1:samp_per_sym) = in(N+1:end, 1);
            
            % Other symbols overlap with previous
            for i = 2:K
                current_index = (i-1)*samp_per_sym + 1 - N;
                out(current_index:current_index+samp_per_sym+N-1) = out(current_index:current_index+samp_per_sym+N-1) + in(:, i);
            end
        end
        
        function out = add_windows(obj, in)
            out = in;
            N = obj.window_length;
            if obj.use_rrc_window
                out(1:N) = in(1:N) .* obj.rrc_taps;
                out(end-N+1:end) = in(end-N+1:end) .* flip(obj.rrc_taps);
            end
        end
        
        function generate_rrc(obj)
            N = obj.window_length;
            obj.rrc_taps = zeros(N, 1);
            for i = 1 : N
                obj.rrc_taps(i) =  0.5 * (1 - sin(pi*(N + 1 - 2*i)/(2*N)));
            end
        end
        
        function out = add_cyclic_prefic(obj, in)
            total_cp_length = obj.cp_length + obj.window_length;
            out = zeros(length(in) + total_cp_length, 1);
            out(total_cp_length + 1:end) = in;
            out(1:total_cp_length) = in(end - total_cp_length+1:end);
        end
        
        
        function out = remove_cyclic_prefix(obj, in)
            out = in(obj.cp_length+1:end);
        end
        
        
        function out = frequency_to_time_domain(obj, in)
            %frequency_to_time_domain Method that can be used to perform the
            %IDFT to go to the time domain signal for OFDM. It will
            %zero pad according to LTE standards. No CP is used.
            %
            % Args:
            %     in:  vector of subcarriers to perform IDFT on.
            %
            %	Author:	Chance Tarver (2018)
            %		tarver.chance@gmail.com
            %
            
            ifft_input = zeros(obj.fft_size, 1);
            ifft_input(2:obj.n_subcarriers/2 + 1) = in(obj.n_subcarriers/2 + 1:end);
            ifft_input(end - obj.n_subcarriers/2 + 1 :end) = in(1:obj.n_subcarriers/2);
            out = ifft(ifft_input);
        end
        
        
        function out = demod(obj, in)
            %demod. Method to demodulate the time domain signal back into the
            %original data.
            %
            % This assumes the input data here is the same length as the
            % original and that everything is synchronized properly.
            
            resource_grid = zeros(obj.n_subcarriers, obj.n_symbols);
            symbol_length = obj.fft_size + obj.cp_length;
            
            for i = 0:obj.n_symbols-1
                td_symbol = in(symbol_length*i + 1: symbol_length*(i+1));
                td = obj.remove_cyclic_prefix(td_symbol);
                fd = obj.time_domain_to_frequency(td);
                resource_grid(:,i+1) = fd;
            end
            out = resource_grid;
        end
        
        
        function evm = calculate_evm(obj, rx_symbols)
            error = rx_symbols - obj.last_tx_freq_d;
            error = error(:); %Convert to vector
            
            original_signal = obj.last_tx_freq_d(:);
            evm = 100 * norm(error) / norm(original_signal);
        end
        
        
        function plot_rx_errors(obj, rx_symbols)
            error = rx_symbols - obj.last_tx_freq_d;
            figure(20);
            hold on;
            title('Error per Subcarrier')
            plot(abs(error));
            xlabel('Subcarrier Index');
            ylabel('abs(error)');
        end
        
        
        function out = time_domain_to_frequency(obj, in)
            %time_domain_to_frequency Method that can be used to perform the
            %DFT to go to the frequency domain signal for OFDM. It assumes
            %zero padding according to LTE standards. No CP is used. Only the
            %occupied subcarriers are returned.
            %
            % Args:
            %     in:  time domain vector to perform DFT on.
            %
            %	Author:	Chance Tarver (2018)
            %		tarver.chance@gmail.com
            %
            
            fftout = fft(in);
            
            out = zeros(obj.n_subcarriers, 1);
            out(1:obj.n_subcarriers/2) = fftout(end - obj.n_subcarriers/2 + 1:end);
            out(obj.n_subcarriers/2+1:end) = fftout(2:obj.n_subcarriers/2+1);
        end
        
        
        function alphabet = QAM_Alphabet(obj, ModulationType)
            %QAM_Alphabet Function to create an alphabet of points of the
            %constellation
            
            %Set up some dictionaries
            modulation_dictionary = containers.Map(obj.constellation_library, ...
                obj.constellaton_order);
            
            MQAM = modulation_dictionary(ModulationType);
            
            alphaMqam = -(sqrt(MQAM)-1):2:(sqrt(MQAM)-1);
            A = repmat(alphaMqam,sqrt(MQAM),1);
            B = flipud(A');
            const_qam = A+1j*B;
            const_qam = const_qam(:);
            alphabet = const_qam;
        end
    end
end
