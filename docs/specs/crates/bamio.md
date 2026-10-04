# bamio — BAM records, compression, and sorting

The crate owns record encoding, sequential BAM inspection, BGZF writing, and coordinate sorting. Candidate discovery and scoring use indexed rust-htslib access separately.

## Record and coordinate contract

`OutRecord` contains a raw-byte read name, SAM flags, MAPQ, optional `AlnView`, sequence, Phred qualities, and optional mate coordinates. `AlnView` contains a contig identifier, 0-based position, CIGAR, and ordered auxiliary tags. CIGAR operations are `Match`, `Ins`, `Del`, `RefSkip`, and `SoftClip`. Tags support character, integer, and string values.

Mapped minus-strand records store reverse-complemented sequence and reversed qualities. Plus-strand and unmapped records retain their original read orientation. Reverse complementation recognizes A/C/G/T case-insensitively, emits uppercase bases, and maps other symbols to N. `restore_original` reverses the mapped-strand transformation when required.

Contig identifiers follow header SQ order. Qualities contain raw Phred values; an empty quality vector represents missing qualities and is encoded as `0xFF`. A nonempty quality vector must match sequence length. Unmapped records use the unmapped flag and have no aligned reference view.

## Encoding and writing

The direct encoder writes standard BAM record fields and falls back to the noodles record representation when the record cannot use the direct encoding path. Sequence uses SAM four-bit encoding, with the first base in the high nibble. Ordered tags are emitted in the supplied order.

`build_header` preserves the supplied contig order and comment. The BGZF writer uses a bounded compression worker count and writes blocks in input order. Record order is supplied by the caller.

## Reading

`open_sequential` returns a contig table and a record iterator. `InRecord` exposes the stored sequence, qualities, flags, reference and mate information available to the parser, CIGAR, and supported auxiliary tags. Stored sequence is interpreted in its BAM orientation. Malformed magic, truncated records, unsupported CIGAR operations, and I/O failures are errors.

## Coordinate sorting

`coordinate_sort(input, output, options)` stably sorts mapper-produced BAM records by contig and position. Unmapped records follow mapped records. Equal-coordinate records retain their input order.

The sorter reads bounded chunks, sorts each chunk, and performs a stable merge. `SortOptions` controls the in-memory record limit, temporary directory, and compression threads. Temporary chunks are removed on success. The output header is preserved, and a `.bai` index is built after writing.

The sorting contract covers the mapper's CIGAR operations and auxiliary tag types. Generic sorting of arbitrary third-party BAM records is outside this interface.

## Validation invariants

- Encoding and decoding preserve supported record fields and tag order.
- Minus-strand sequence and qualities follow the SAM orientation rule.
- Sorting preserves the record multiset and yields nondecreasing coordinates.
- Changing chunk boundaries does not alter record order or content.
- The generated BAM index supports region access.
