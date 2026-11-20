# PTAPI version 1

The core of version 1 (`include/`) and the problem used to measure it (`SD/`,
with one instance), unchanged. They are here only so that `bench/` can compare
the speed of the two versions on the same machine: `make bench` and
`make bench-sd` in `../bench`.

New programs should use the headers of `../include` (version 2). Code written
for version 1 compiles with them without changes; see
[docs/migrating-from-v1.md](../docs/migrating-from-v1.md).
