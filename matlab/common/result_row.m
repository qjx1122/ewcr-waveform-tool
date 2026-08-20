function row = result_row(r,dataset,signal,param)
%RESULT_ROW Convert unified result struct to one scalar row struct.
%
% v2.4 compatibility rule:
%   - text fields are stored as plain char vectors in each scalar struct;
%   - numeric fields are scalar numeric values;
%   - conversion of a struct array to a MATLAB table is handled ONLY by
%     common/structrows_to_table.m.
%
% This deliberately avoids both string() and struct2table()'s implicit
% cell/char expansion rules.

if ~isstruct(r) || ~isscalar(r)
    error('FiveCodecBenchmark:NonScalarUnifiedResult', ...
        ['result_row expected one scalar unified codec result, but received ' ...
         '%d struct elements. Check run_codec_unified construction.'], numel(r));
end

m = r.metrics;
z = r.resource;

row = struct();

% Text metadata as plain character vectors.
row.dataset = char(dataset);
row.signal  = char(signal);
row.algo    = char(r.algo);
row.param   = char(param);

% Signal / rate-distortion metrics.
row.N = m.N;
row.CR = m.CR;
row.rate_bps = m.rate_bps;
row.NMSE_dB = m.NMSE_dB;
row.SNR_dB = m.SNR_dB;
row.RMSE = m.RMSE;
row.PRD = m.PRD;
row.MAXE = m.MAXE;
row.PSNR_dB = m.PSNR_dB;
row.fund_amp_relerr = m.fund_amp_relerr;
row.fund_phase_err_deg = m.fund_phase_err_deg;
row.THD_abs_error = m.THD_abs_error;

% Measured codec time for this invocation.
row.enc_time_s = r.enc_time_s;
row.dec_time_s = r.dec_time_s;

% Resource / payload accounting.
row.logical_payload_bytes = z.logical_payload_bytes;
row.matlab_packet_bytes = z.matlab_packet_bytes;
row.encoder_working_est_bytes = z.encoder_working_est_bytes;
row.decoder_working_est_bytes = z.decoder_working_est_bytes;
row.edge_state_est_bytes = z.edge_state_est_bytes;
row.source_mfile_bytes = z.source_mfile_bytes;
end
