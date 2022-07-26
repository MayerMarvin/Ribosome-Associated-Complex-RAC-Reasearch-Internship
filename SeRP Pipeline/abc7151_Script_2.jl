#!/usr/bin/env julia

Base.CoreLogging.disable_logging(Base.CoreLogging.Debug)

using BioAlignments
using BioSequences
using GenomicFeatures
using DataStructures

using CodecZlib
using HDF5

using Gadfly
import Cairo
import Compose

using ArgParse

@enum TranscriptChoice longestCDS longestTranscript
@enum AssignmentType FivePrime ThreePrime
const ProteinAttributes = Dict{String, String}("protein_id" => "protein_id",
                                                "Name" => "protein_name",
                                                "orf_classification" => "orf_classification")
const mRNAAttributes = Dict{String, String}("transcript_id" => "transcript_id",
                                            "Name" => "transcript_name")
const GeneAttributes = Dict{String, String}("gene" => "gene")

struct ReadLengthRange
    start::Union{Unsigned, Nothing}
    stop::Union{Unsigned, Nothing}
end
ReadLengthRange() = ReadLengthRange(nothing, nothing)
function ReadLengthRange(x::AbstractString)
    ranges = split(x, ':', keepempty=true)
    if length(ranges) > 2
        throw(ArgumentError("\"$x\" is not a valid read length range."))
    end
    start = length(ranges[1]) > 0 ? parse(UInt, ranges[1]) : Nothing
    stop = length(ranges[2]) > 0 ? parse(UInt, ranges[2]) : Nothing
    ReadLengthRange(start, stop)
end

function assign3prime(read::Union{SAM.Record, BAM.Record}, clip5::Unsigned)
    fmt = parentmodule(typeof(read))
    aln = fmt.alignment(read)
    ispositivestrand = fmt.flag(read) & 0x10 == 0
    seqlen = fmt.seqlength(read)
    offsets = [13, 16, 19]
    if seqlen < maximum(offsets)
        return nothing
    end

    assignments = [UInt(k) for (k, nothing) in (ispositivestrand ? seq2ref.((aln,), seqlen .- offsets .+ 1) : seq2ref.((aln,), offsets))]
    return assignments
end

function assign5prime(read::Union{SAM.Record, BAM.Record}, clip5::Unsigned)
    fmt = parentmodule(typeof(read))
    aln = fmt.alignment(read)
    ispositivestrand = fmt.flag(read) & 0x10 == 0
    seqlen = fmt.seqlength(read)
    offsets = [16, 13, 10]
    if clip5 > 0
        if ispositivestrand
            rng = 1:clip5
        else
            rng = seqlen:seqlen - clip5 + 1
        end
        for i ∈ rng
            if seq2ref(aln, i)[2] == OP_SOFT_CLIP
                offsets .+= 1
            else
                break
            end
        end
    end
    if seqlen < maximum(offsets)
        return nothing
    end

    assignments = [UInt(k) for (k, nothing) in (ispositivestrand ? seq2ref.((aln,), offsets) : seq2ref.((aln,), seqlen .-offsets .+ 1))]
    return assignments
end

@enum RibosomeSite asite psite esite

abstract type AbstractRead end

struct Read <: AbstractRead
    refname::String
    strand::Strand
    readlength::UInt16
    clippedlength::UInt16
    lmap::UInt32
    rmap::UInt32

    asite::UInt32
    psite::UInt32
    esite::UInt32
end

refname(rd::Read) = rd.refname
strand(rd::Read) = rd.strand
rdlength(rd::Read) = rd.readlength
clippedlength(rd::Read) = rd.clippedlength
lmap(rd::Read) = rd.lmap
rmap(rd::Read) = rd.rmap
function site(rd::Read, s::RibosomeSite)
    if s == asite
        rd.asite
    elseif s == psite
        rd.psite
    elseif s == esite
        rd.esite
    end
end

struct UMIRead <: AbstractRead
    umi::DNAKmer
    read::Read
end

umiseq(rd::UMIRead) = rd.umi
refname(rd::UMIRead) = refname(rd.read)
strand(rd::UMIRead) = strand(rd.read)
rdlength(rd::UMIRead) = rdlength(rd.read)
clippedlength(rd::UMIRead) = clippedlength(rd.read)
lmap(rd::UMIRead) = lmap(rd.read)
rmap(rd::UMIRead) = rmap(rd.read)
site(rd::UMIRead, s::RibosomeSite) = site(rd.read, s)

function Read(record::Union{SAM.Record, BAM.Record}, clip5::Unsigned, asstype::AssignmentType)
    fmt = parentmodule(typeof(record))
    aln = fmt.alignment(record)
    seq = fmt.sequence(record)
    seqlen = fmt.seqlength(record)
    clippedlen = seqlen
    for i ∈ 1:clip5
        if seq2ref(aln, i)[2] == OP_SOFT_CLIP
            clippedlen -= 1
        else
            break
        end
    end
    ispositivestrand = fmt.flag(record) & 0x10 == 0
    ass = asstype == FivePrime ? assign5prime(record, clip5) : assign3prime(record, clip5)
    if ass === nothing
        return nothing
    else
        return Read(fmt.refname(record),
                    ispositivestrand ? STRAND_POS : STRAND_NEG,
                    seqlen,
                    clippedlen,
                    seq2ref(aln, 1)[1],
                    seq2ref(aln, seqlen)[1],
                    ass...)
    end
end

function UMIRead(record::Union{SAM.Record, BAM.Record}, clip5::Unsigned, asstype::AssignmentType)
    fmt = parentmodule(typeof(record))
    umistr = fmt.tempname(record)
    umiseq = DNASequence(umistr[(findlast("____", umistr).stop + 1):end])

    aln = fmt.alignment(record)
    seq = fmt.sequence(record)
    seqlen = fmt.seqlength(record)
    ispositivestrand = fmt.flag(record) & 0x10 == 0
    for i ∈ 1:clip5
        if seq2ref(aln, i)[2] == OP_SOFT_CLIP
            push!(umiseq, seq[i])
        else
            break
        end
    end
    if findfirst(DNA_N, umiseq) != nothing
        return nothing
    end
    rd = Read(record, clip5, asstype)
    if rd === nothing
        return nothing
    else
        return UMIRead(DNAKmer(umiseq), rd)
    end
end

function Base.hash(x::UMIRead, h::UInt)
    hash(umiseq(x),
         hash(refname(x),
              hash(strand(x),
                   hash(rdlength(x),
                        hash(lmap(x),
                             hash(rmap(x), h))))))
end

function Base.:(==)(x::UMIRead, y::UMIRead)
    umiseq(x) == umiseq(y) &&
    refname(x) == refname(y) &&
    strand(x) == strand(y) &&
    rdlength(x) == rdlength(y) &&
    lmap(x) == lmap(y) &&
    rmap(x) == rmap(y)
end

const UMICounts = DefaultDict{String, DefaultDict{UInt32, UInt32}}
UMICounts() = UMICounts(() -> DefaultDict{UInt32, UInt32}(UInt32(0)))

function featurelength(ic::IntervalCollection{T}) where T
    len = UInt64(0)
    for i in ic
        len += rightposition(i) - leftposition(i) + 1
    end
    len
end

struct mRNA
    parent::String
    strand::GenomicFeatures.Strand
    ID::String
    chromosome::String
    transcript::Interval{GFF3.Record}
    attributes::Dict{String, Union{String, Number}}
end

mRNA(parent::String, strand::GenomicFeatures.Strand, ID::String, chromosome::String, transcript::Interval{GFF3.Record}) = mRNA(parent, strand, ID, chromosome, transcript, Dict{String, Union{String, Number}}())
mRNA(parent::String, strand::GenomicFeatures.Strand, ID::String, chromosome::String) = mRNA(parent, strand, ID, chromosome, Interval{GFF3.Record}())
Base.length(p::mRNA) = rightposition(p.transcript) - leftposition(p.transcript) + 1

struct Protein
    parent::String
    strand::GenomicFeatures.Strand
    ID::String
    chromosome::String
    CDS::IntervalCollection{GFF3.Record}
    attributes::Dict{String, Union{String, Number}}
end

Protein(transcript::String, strand::GenomicFeatures.Strand, chromosome::String, ID::String) = Protein(transcript, strand, chromosome, ID, IntervalCollection{GFF3.Record}(), Dict{String, Union{String, Number}}())
Base.length(p::Protein) = featurelength(p.CDS)
Base.push!(p::Protein, i::Interval{T}) where T = push!(p.CDS, i)

function Base.iterate(p::Protein)
    if length(p.CDS) == 0
        nothing
    else
        iterate(p, (iterate(p.CDS), (0, p.strand == STRAND_POS ? UInt64(1) : length(p))))
    end
end
function Base.iterate(p::Protein, state)
    (int, intState), (pos, cdspos) = state
    genomicPos = leftposition(int) + pos
    if genomicPos > rightposition(int)
        pos = 0
        nstate = iterate(p.CDS, intState)
        if nstate === nothing
            return nothing
        else
            int, intState = nstate
            genomicPos = leftposition(int) + pos
        end
    end
    ((genomicPos, cdspos), ((int, intState), (pos + 1, p.strand == STRAND_POS ? cdspos + 1 : cdspos - 1)))
end

struct Gene
    name::String
    ID::String
    chromosome::String
    mRNAs::Dict{String, mRNA}
    proteins::Dict{String, Protein}
    attributes::Dict{String, Union{String, Number}}
end

Gene(name::String, ID::String, chromosome::String) = Gene(name, ID, chromosome, Dict{String, mRNA}(), Dict{String, Protein}(), Dict{String, Union{String, Number}}())

const StrandPair = NamedTuple{(:fw, :rev), Tuple{DefaultDict{UInt32, UInt32, UInt32}, DefaultDict{UInt32, UInt32, UInt32}}}
StrandPair() = StrandPair((DefaultDict{UInt32, UInt32, UInt32}(UInt32(0)), DefaultDict{UInt32, UInt32, UInt32}(UInt32(0))))
strandcounts(x::StrandPair, s::Strand) = s == STRAND_POS ? x[:fw] : x[:rev]

@enum SiteCountsType all umi
abstract type AbstractSiteCounts end

struct UMISiteCounts <: AbstractSiteCounts
    all::StrandPair
    umi::StrandPair
end
UMISiteCounts() = UMISiteCounts(StrandPair(), StrandPair())
types(::Type{UMISiteCounts}) = [all, umi]
types(x::UMISiteCounts) = types(UMISiteCounts)

struct SiteCounts <: AbstractSiteCounts
    all::StrandPair
end
SiteCounts() = SiteCounts(StrandPair())
types(::Type{SiteCounts}) = [all]
types(x::SiteCounts) = types(SiteCounts)

counts(sitecounts::AbstractSiteCounts, type::SiteCountsType)::StrandPair = getfield(sitecounts, Symbol(type))

const GenomicPositions{T} = DefaultDict{String, T} where T <: AbstractSiteCounts
GenomicPositions(umis::Bool) = umis ? GenomicPositions{UMISiteCounts}(() -> UMISiteCounts()) : GenomicPositions{SiteCounts}(() -> SiteCounts())
types(p::GenomicPositions) = types(valtype(p))

struct GenomicSitePositions{T}
    asite::GenomicPositions{T}
    psite::GenomicPositions{T}
    esite::GenomicPositions{T}
end
GenomicSitePositions(umis::Bool) = GenomicSitePositions(GenomicPositions(umis), GenomicPositions(umis), GenomicPositions(umis))
function site(p::GenomicSitePositions, s::RibosomeSite)
    if s == asite
        p.asite
    elseif s == psite
        p.psite
    elseif s == esite
        p.esite
    end
end
function types(p::GenomicSitePositions)
    tps = types(p.asite)
    @assert types(p.psite) == types(p.esite) == tps
    tps
end

abstract type AbstractFootprintLengthCounts end

struct UMIFootprintLengthCounts <: AbstractFootprintLengthCounts
    all::DefaultDict{UInt8, UInt32}
    accepted::DefaultDict{UInt8, UInt32}
    acceptedumi::DefaultDict{UInt8, UInt32}
end
UMIFootprintLengthCounts() = UMIFootprintLengthCounts(DefaultDict{UInt8, UInt32}(UInt32(0)), DefaultDict{UInt8, UInt32}(UInt32(0)), DefaultDict{UInt8, UInt32}(UInt32(0)))

struct FootprintLengthCounts <: AbstractFootprintLengthCounts
    all::DefaultDict{UInt8, UInt32}
    accepted::DefaultDict{UInt8, UInt32}
end
FootprintLengthCounts() = FootprintLengthCounts(DefaultDict{UInt8, UInt32}(UInt32(0)), DefaultDict{UInt8, UInt32}(UInt32(0)))

function testGZIP(in::IOStream)
    pos = position(in)
    if eof(in)
        seekstart(in)
    end

    magicnr = read(in, 2)
    if length(magicnr) < 2
        return false
    end
    seek(in, pos)
    return magicnr[1] == 0x1f && magicnr[2] == 0x8b
end

function xread(infile::AbstractString)
    in = open(infile, "r")
    return testGZIP(in) ? GzipDecompressorStream(in) : in
end

function parseGFF3(infile::AbstractString)
    @info "reading genome annotation"
    proteins = Dict{String, Protein}()
    mRNAs = Dict{String, mRNA}()
    genes = Dict{String, Gene}()
    reader = GFF3.Reader(xread(infile))
    for record in reader
        if GFF3.featuretype(record) == "CDS"
            att = Dict(GFF3.attributes(record))
            id = haskey(att, "ID") ? att["ID"][1] : att["Name"][1]
            chrom = GFF3.seqid(record)
            if haskey(proteins, id)
                prot = proteins[id]
            else
                prot = Protein(att["Parent"][1], GFF3.strand(record), id, chrom)
                for k ∈ keys(ProteinAttributes)
                    if haskey(att, k)
                        prot.attributes[k] = att[k][1]
                    end
                end
                proteins[id] = prot
            end
            push!(prot, convert(Interval, record))
        elseif GFF3.featuretype(record) == "mRNA"
            att = Dict(GFF3.attributes(record))
            id = haskey(att, "ID") ? att["ID"][1] : att["Name"][1]
            chrom = GFF3.seqid(record)
            if haskey(mRNAs, id)
                mrna = mRNAs[id]
            else
                mrna = mRNA(att["Parent"][1], GFF3.strand(record), id, chrom, convert(Interval, record))
                for k ∈ keys(mRNAAttributes)
                    if haskey(att, k)
                        mrna.attributes[k] = att[k][1]
                    end
                end
                mRNAs[id] = mrna
            end
        elseif GFF3.featuretype(record) == "gene"
            att = Dict(GFF3.attributes(record))
            id = haskey(att, "ID") ? att["ID"][1] : att["Name"][1]
            chrom = GFF3.seqid(record)
            gene = Gene(att["Name"][1], id, chrom)
            for k ∈ keys(GeneAttributes)
                    if haskey(att, k)
                        gene.attributes[k] = att[k][1]
                    end
                end
            genes[id] = gene
        end
    end
    close(reader)
    for (id, mrna) ∈ mRNAs
        if haskey(genes, mrna.parent)
            genes[mrna.parent].mRNAs[id] = mrna
        end
    end
    for (id, prot) ∈ proteins
        if haskey(mRNAs, prot.parent)
            genes[mRNAs[prot.parent].parent].proteins[id] = prot
        elseif haskey(genes, prot.parent)
            genes[prot.parent].proteins[id] = prot
        end
    end
    for k ∈ [id for (id, gene) ∈ genes if isempty(gene.mRNAs) && isempty(gene.proteins)]
        delete!(genes, k)
    end

    return genes
end

function processUMIGroup(calignments::Accumulator{UMIRead, <:Unsigned}, positions::GenomicSitePositions, ucounts::UMICounts, lcounts::UMIFootprintLengthCounts)
    for (umird, count) ∈ calignments
        ucounts[refname(umird)][count] += 1
        lcounts.acceptedumi[clippedlength(umird)] += 1
        lcounts.accepted[clippedlength(umird)] += count
        for rsite ∈ instances(RibosomeSite)
            for (type, cnt) ∈ zip((all, umi), (count, 1))
                strandcounts(counts(site(positions, rsite)[refname(umird)], type), strand(umird))[site(umird, rsite)] += cnt
            end
        end
    end
    empty!(calignments.map)
    nothing
end

function readUMIAlignments(reader::Union{BAM.Reader, SAM.Reader}, clip5::Unsigned=0, lengthrange::ReadLengthRange=ReadLengthRange(), asstype::AssignmentType=FivePrime)
    hd = findall(header(reader), "HD")[1]
    hdvals = Dict(SAM.keyvalues(hd))
    if !haskey(hdvals, "SO") || hdvals["SO"] != "coordinate"
        @error "alignment file must be coordinate-sorted"
    end

    positions = GenomicSitePositions(true)
    counts = UMICounts()
    lengthcounts = UMIFootprintLengthCounts()
    calignments = Accumulator{UMIRead, UInt32}()
    cpos = UInt32(0)
    cref = ""
    record = parentmodule(typeof(reader)).Record()
    while !eof(reader)
        read!(reader, record)
        rd = UMIRead(record, clip5, asstype)
        if rd != nothing
            lengthcounts.all[clippedlength(rd)] += 1
            if lengthrange.start !== nothing && rdlength(rd) < lengthrange.start || lengthrange.stop !== nothing && rdlength(rd) > lengthrange.stop
                continue
            elseif lmap(rd) == cpos && refname(rd) == cref
                push!(calignments, rd)
            else
                processUMIGroup(calignments, positions, counts, lengthcounts)
                cref = refname(rd)
                cpos = lmap(rd)
                push!(calignments, rd)
            end
        end
    end
    processUMIGroup(calignments, positions, counts, lengthcounts)
    return positions, lengthcounts, counts
end

function readAlignments(reader::Union{BAM.Reader, SAM.Reader}, clip5::Unsigned=0, lengthrange::ReadLengthRange=ReadLengthRange(), asstype::AssignmentType=FivePrime)
    positions = GenomicSitePositions(false)
    lengthcounts = FootprintLengthCounts()
    record = parentmodule(typeof(reader)).Record()
    while !eof(reader)
        read!(reader, record)
        rd = Read(record, clip5, asstype)
        if rd != nothing
            lengthcounts.all[clippedlength(rd)] += 1
            if lengthrange.start !== nothing && rdlength(rd) < lengthrange.start || lengthrange.stop !== nothing && rdlength(rd) > lengthrange.stop
                continue
            end
            lengthcounts.accepted[rdlength(rd)] += 1
            for rsite ∈ instances(RibosomeSite)
                strandcounts(counts(site(positions, rsite)[refname(rd)], all), strand(rd))[site(rd, rsite)] += 1
            end
        end
    end
    return positions, lengthcounts, nothing
end

function readAlignments(infile::AbstractString, clip5::Unsigned, lengthrange::ReadLengthRange=ReadLengthRange(), use_umis::Bool=false, asstype::AssignmentType=FivePrime)
    @info "reading alignments"
    in = open(infile, "r")
    fmt = testGZIP(in) ? BAM : SAM
    reader = fmt.Reader(in)::Union{BAM.Reader, SAM.Reader}

    if use_umis
        positions, lengthcounts, counts = readUMIAlignments(reader, clip5, lengthrange, asstype)
    else
        positions, lengthcounts, counts = readAlignments(reader, clip5, lengthrange, asstype)
    end
    close(reader)
    return positions, lengthcounts, counts
end

function writeGenomeCoordinates(positions::GenomicSitePositions, outdir::AbstractString)
    out_counts = Dict((i, h5open(joinpath(outdir, "counts_" * string(i) * ".h5"), "w")) for i ∈ types(positions))
    for rsite ∈ instances(RibosomeSite)
        position = site(positions, rsite)
        g = Dict((i, g_create(f, string(rsite))) for (i, f) ∈ out_counts)
        for (chrom, scounts) ∈ position
            for (type, h5g) ∈ g
                counts_fw = strandcounts(counts(scounts, type), STRAND_POS)
                counts_rev = strandcounts(counts(scounts, type), STRAND_NEG)

                occupiedpos = sort(collect(union(keys(counts_fw), keys(counts_rev))))
                pos, fstrand, rstrand = Vector{UInt32}(),
                                        Vector{UInt32}(),
                                        Vector{UInt32}()
                for p ∈ occupiedpos
                    push!(pos, p)
                    push!(fstrand, p ∈ keys(counts_fw) ? counts_fw[p] : 0) # DefaultDict inserts the default value into the dict already upon retrieval,
                    push!(rstrand, p ∈ keys(counts_rev) ? counts_rev[p] : 0) # this prevents filling the Dicts up with zeroes
                end
                h5g[chrom, "compress", 9] = hcat(pos, fstrand, rstrand)
            end
        end
    end
    for f ∈ values(out_counts)
        close(f)
    end
end

function chooseLongestCDS(CDSs, genetrx)
    lengths = [(length(prot), haskey(genetrx, prot.parent) ? length(genetrx[prot.parent]) : 0, id) for (id, prot) ∈ CDSs]
    sort!(lengths, rev=true)
    [CDSs[id] for (nothing, nothing, id) ∈ lengths]
end

function chooseLongestTranscript(CDSs, genetrx)
    lengths = [(haskey(genetrx, prot.parent) ? length(genetrx[prot.parent]) : 0, length(prot), id) for (id, prot) ∈ CDSs]
    sort!(lengths, rev=true)
    [CDSs[id] for (nothing, nothing, id) ∈ lengths]
end

function aggregateAlignments(positions::GenomicSitePositions, genes::AbstractDict{<:AbstractString, Gene}, outdir::AbstractString; trxChoice = longestCDS)
    @info "assigning counts to genes"
    for rsite ∈ instances(RibosomeSite)
        position = site(positions, rsite)
        ptypes = types(position)
        out = Dict((type, h5open(joinpath(outdir, "cds_" * string(type) * "_" * string(rsite) * ".h5"), "w")) for type ∈ ptypes)
        pos = Vector{UInt32}()
        cds_counts = Dict((type, Vector{UInt32}()) for type ∈ ptypes)
        for gene ∈ values(genes)
            if trxChoice == longestCDS
                proteins = chooseLongestCDS(gene.proteins, gene.mRNAs)
            elseif trxChoice == longestTranscript
                proteins = chooseLongestTranscript(gene.proteins, gene.mRNAs)
            end

            for protein ∈ proteins
                if haskey(position, protein.chromosome)
                    gcounts = position[protein.chromosome]
                    for (genomicCoord, cdsCoord) ∈ protein
                        have_counts = false
                        for type ∈ ptypes
                            scounts = strandcounts(counts(gcounts, type), protein.strand)
                            if genomicCoord ∈ keys(scounts)
                                push!(cds_counts[type], scounts[genomicCoord])
                                have_counts = true
                            end
                        end
                        if have_counts
                            push!(pos, cdsCoord)
                        end
                    end
                    if length(pos) > 0
                        for (type, h5out) ∈ out
                            arr = hcat(pos, cds_counts[type])
                            dset = d_create(h5out, gene.name, arr)[1]
                            dset[:,:] = arr
                            attrs(dset)["chromosome"] = protein.chromosome
                            attrs(dset)["cds_length"] = length(protein)
                            attrs(dset)["strand"] = string(protein.strand)

                            for (k, a) ∈ protein.attributes
                                attrs(dset)[ProteinAttributes[k]] = a
                            end
                            if  haskey(gene.mRNAs, protein.parent)
                                for (k, a) ∈ gene.mRNAs[protein.parent].attributes
                                    attrs(dset)[mRNAAttributes[k]] = a
                                end
                            end
                            for (k, a) ∈ gene.attributes
                                attrs(dset)[GeneAttributes[k]] = a
                            end
                        end
                    end
                    empty!(pos)
                    for cnt ∈ values(cds_counts)
                        empty!(cnt)
                    end
                    break
                end
            end
        end
        for f ∈ values(out)
            close(f)
        end
    end
    nothing
end

function aggregateAlignments(positions::GenomicSitePositions, gfffile::AbstractString, outdir::AbstractString; trxChoice = longestCDS)
    genes = parseGFF3(gfffile)
    aggregateAlignments(positions, genes, outdir, trxChoice=trxChoice)
end

function plotLengthCounts(counts::AbstractDict{<:Unsigned, <:Unsigned}, out::Cairo.CairoSurface, title::Union{Nothing, AbstractString})
    c = Cairo.CairoContext(out)
    s = Compose.CAIROSURFACE(out, c)
    xrange = extrema(keys(counts))
    Compose.draw(s, plot(x=collect(keys(counts)),
                         y=collect(values(counts)),
                         Geom.bar,
                         Guide.title(title),
                         Guide.xlabel("read length / nucleotides"),
                         Guide.ylabel("count"),
                         Scale.x_continuous(minticks=4),
                         Coord.cartesian(xmin=xrange[1], xmax=xrange[2])))
    Cairo.show_page(c)
end

function plotLengthCounts(counts::UMIFootprintLengthCounts, outdir::AbstractString)
    out =  Cairo.CairoPDFSurface(joinpath(outdir, "footprint_length_distribution.pdf"), 5*72, 5*72)
    plotLengthCounts(counts.all, out, "all reads")
    plotLengthCounts(counts.accepted, out, "accepted reads")
    plotLengthCounts(counts.acceptedumi, out, "accepted, UMI-collapsed reads")
    Cairo.finish(out)
    nothing
end

function plotLengthCounts(counts::FootprintLengthCounts, outdir::AbstractString)
    out =  Cairo.CairoPDFSurface(joinpath(outdir, "footprint_length_distribution.pdf"), 5*72, 5*72)
    plotLengthCounts(counts.all, out, "all reads")
    plotLengthCounts(counts.accepted, out, "accepted reads")
    Cairo.finish(out)
    nothing
end

function plotUMICounts(counts::UMICounts, outdir::AbstractString)
    @info "plotting reads/UMI distributions"

    out =  Cairo.CairoPDFSurface(joinpath(outdir, "reads_per_umi.pdf"), 5*72, 5*72)
    for (chromosome, cnt) ∈ counts
        c = Cairo.CairoContext(out)
        s = Compose.CAIROSURFACE(out, c)
        Compose.draw(s, plot(x=collect(keys(cnt)),
                             y=collect(values(cnt)),
                             Geom.bar,
                             Scale.y_log10,
                             Guide.title(chromosome),
                             Guide.xlabel("reads / UMI"),
                             Guide.ylabel("count")))
        Cairo.show_page(c)
    end
    Cairo.finish(out)
    nothing
end

function Base.show(io::IO, rng::ReadLengthRange)
    print(io, rng.start === nothing ? "" : Int(rng.start))
    print(io, ':')
    print(io, rng.stop === nothing ? "" : Int(rng.stop))
end

function ArgParse.parse_item(::Type{TranscriptChoice}, x::AbstractString)
    validOptions = instances(TranscriptChoice)
    opt = findfirst(string.(validOptions) .== x)
    if opt === nothing
        throw(ArgumentError("invalid transcript choice \"" * x * "\""))
    else
        validOptions[opt]
    end
end

function ArgParse.parse_item(::Type{AssignmentType}, x::AbstractString)
    if x == "5"
        FivePrime
    elseif x == "3"
        ThreePrime
    else
        throw(ArgumentError("invalid assignment type \"" * x * "\""))
    end
end

function ArgParse.parse_item(::Type{ReadLengthRange}, x::AbstractString)
    ReadLengthRange(x)
end

s = ArgParseSettings()
@add_arg_table s begin
    "--gff", "-g"
        help = "Path to genome annotation in GFF3 format (may be GZIP-compressed)."
        arg_type = String
        required = true
    "--clip5", "-c"
        help = "Number of bases at the 5' end that may be generated by untemplated addition."
        arg_type = UInt8
        default = UInt8(0)
        required = false
    "--use_umis", "-u"
        help = "Use UMI information to discard PCR duplicates. This option assumes that the UMI is stored as the last part of the read name, separated by '____' from the actual read ID. INFILE must be coordinate-sorted."
        action = :store_true
    "--length_range", "-l"
        help = "Valid read lengths. Must be a range expression, e.g. \"20:50\" will only include reads that are between 20 and 50 bases long (after 5' clipping), whereas \"20:\" will include all reads that are at least 20 bases long."
        arg_type = ReadLengthRange
        default = ReadLengthRange()
        required = false
    "--transcript_choice", "-t"
        help = "How to choose the reported transcript for CDS assignments. Valid choices are " * join(string.(instances(TranscriptChoice)), " ,", ", or ")
        arg_type = TranscriptChoice
        default = longestCDS
        required = false
    "--assignment_type", "-a"
        help = "Ribosome site assignment type. Valid choices are 5 for 5'-assignnment or 3 for 3'-assignment."
        arg_type = AssignmentType
        default = FivePrime
        required = false
    "--outdir", "-o"
        help = "output directory"
        arg_type = String
        required = true
    "INFILE"
        help = "path to alignment file in SAM or BAM format"
        arg_type = String
        required = true
end
args = parse_args(s, as_symbols=true)

mkpath(args[:outdir])

positions, lengthCounts, umiCounts = readAlignments(args[:INFILE], args[:clip5], args[:length_range], args[:use_umis], args[:assignment_type])

plotLengthCounts(lengthCounts, args[:outdir])

if args[:use_umis]
    plotUMICounts(umiCounts, args[:outdir])
end

writeGenomeCoordinates(positions, args[:outdir])
aggregateAlignments(positions, args[:gff], args[:outdir], trxChoice=args[:transcript_choice])
