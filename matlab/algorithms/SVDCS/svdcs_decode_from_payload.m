function xhat = svdcs_decode_from_payload(packet)
%SVDCS_DECODE_FROM_PAYLOAD Verify decoding after actual LZW payload recovery.
%
% This function intentionally discards packet.mc and rebuilds MC from the
% compressed byte payload, then invokes the normal spectral decoder.

if packet.cfg.use_lzw
    raw=lzw_decode12(packet.payload,packet.lzw_ncodes);
else
    raw=packet.payload;
end

mc=svdcs_deserialize_mc(raw,2*packet.Nppc,packet.Nframes);
p2=packet;
p2.mc=mc;
xhat=svdcs_decode(p2);
end
