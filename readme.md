# IBSAMRAI2

Backwards-compatible fork of SAMRAI 2.4.4 containing some additional patches,
both for improving performance and for fixing problems with modern compilers.

IBSAMRAI2 requires a C++17 compiler. If the compiler's default language standard
is older, as for GCC 9 and 10, `configure` adds `-std=c++17` to `CXX`; it stops
if the compiler cannot compile C++17. See `./configure --help` for more
information on setting up compilation flags.
