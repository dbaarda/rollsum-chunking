# rollsum-chunking
Experiments in Content Defined Chunking using rollsums.

This was research on chunking algorithms inspired by the FastCDC paper.

See http://minkirri.apana.org.au/wiki/RollsumChunking for early analysis and
observations.

See [RESULTS.rst](RESULTS.rst) for the details and conclusion.

TLDR; FastCDC is not faster and doesn't get better deduplication than the
standard simple exponential chunker for the same average, min, and max chunk
size. Fancier chunk-normalizing algorithms also are not better. The best
standard exponential chunker settings for speed and deduplication are very
different from the generally recommended values.
