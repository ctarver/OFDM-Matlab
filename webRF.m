classdef webRF < handle
    %webRF Class wrapper for the webRF PA.
    % http://dpdcompetition.com/rfweblab/
    
    properties
        RMSin
        RMSout
        Idc
        Vdc
        synchronization
        Fs
    end
    
    methods
        function obj = webRF(dbm_power)
            %webRF Construct an instance of this class
            if nargin == 0
                dbm_power = -24;
            end
            obj.RMSin = dbm_power;
            obj.synchronization.sub_sample = 1;
            obj.Fs = 200e6;
        end
        
        function [y, up, in_up] = transmit(obj, x, current_sample_rate)
            %transmit. Take input signal, x, and broadcast it through the
            %RFWebLab PA.
            %
            %Args:
            %   -x: column vector. Will be normalized in RFWebLab function
            %
            %Returns:
            %   -y: column vector result from sending x through the PA. Y
            %   is normalized to be the same ||.||2 norm as x.
            
            % if there are 2 inputs, we will do up and downsampling too.
            if nargin == 3
                 x = obj.upsample(x, current_sample_rate);
            end
            in_up = x;
            
            if length(x) > 1000000
                warning("Too long for webRF.");
            end
            [y, obj.RMSout, obj.Idc, obj.Vdc] = RFWebLab_PA_meas_v1_1(x, obj.RMSin);
            
            % Need something to guarantee same as input length and aligned in TD.
            y = [y(7:end)];
            length_input = length(x);
            length_output = length(y);
            y = [y; zeros(length_input - length_output, 1)];
            
            % Normalize
            y = y * norm(x) / norm(y);
            if  obj.synchronization.sub_sample
                %Set up a LS estimation for figuring out a subsample delay.
                X = [y [0; y(1:end-1)]];
                coeffs = (X'*X) \ (X'*x);
                y = X*coeffs;
            end
            up = y;
            if nargin == 3
               y = obj.downsample(y, current_sample_rate); 
            end
        end
        
        function out = upsample(obj, in, current_sample_rate)
            desired_sample_rate = obj.Fs;
            upsample_rate = floor(desired_sample_rate/current_sample_rate);
            up = upsample(in, upsample_rate);
            up = [up; zeros(255,1)];
            b = firls(255,[0 (1/upsample_rate -0.02) (1/upsample_rate +0.02) 1],[1 1 0 0]);
            out = filter(b,1,up);
            out = out(128:end);
            out = out(1:length(up)-255);
            out = out * rms(in)/rms(out);
            
        end
        
        function out = downsample(obj, in, desired_sample_rate)
            current_sample_rate = obj.Fs;
            downsample_rate = floor(current_sample_rate/desired_sample_rate);
            out = downsample(in, downsample_rate);
        end
    end
end