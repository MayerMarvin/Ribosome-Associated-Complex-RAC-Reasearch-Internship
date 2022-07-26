#!/usr/bin/env julia

using ArgParse
using CodecZlib

function trimUMIs(in::IO, out::IO, umi3::Unsigned, umi5::Unsigned)
    rd = Vector{String}(undef, 4)
    while !eof(in)
        for i in 1:4
            rd[i] = readline(in)
        end
        if length(rd[2]) <= umi5 + umi3
            continue
        end
        umi = rd[2][1:umi5] * rd[2][end - umi3 + 1:end]

        cutoff = findfirst(isequal(' '), rd[1])
        rd[1] = rd[1][1:(cutoff - 1)] * "____" * umi
        for i in (2,4)
            rd[i] = rd[i][umi5 + 1:end - umi3]
        end
        for l in rd
            write(out, l, '\n')
        end
    end
end

s = ArgParseSettings()
@add_arg_table s begin
    "--umi3"
        help = "number of bases at the 3'-end to use as UMI"
        arg_type = UInt
        default = UInt(0)
        required = false
    "--umi5"
        help = "number of bases at the 5'-end to use as UMI"
        arg_type = UInt
        default = UInt(0)
        required = false
    "INFILE"
        help = "path to FASTQ file (may be GZIP compressed)"
        arg_type = String
        required = true
    "OUTFILE"
        help = "path to output FASTQ file (will be GZIP compressed)"
        arg_type = String
        required = true
end

args = parse_args(s, as_symbols=true)

in = open(args[:INFILE], "r")
magicnr = read(in, 2)
seek(in, 0)
if magicnr[1] == 0x1f && magicnr[2] == 0x8b
    in = GzipDecompressorStream(in)
end

out = GzipCompressorStream(open(args[:OUTFILE], "w"), level=6)
trimUMIs(in, out, args[:umi3], args[:umi5])
close(in)
close(out)
