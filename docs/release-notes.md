# Writing an ibsamrai2 release

1. Recompile in debug mode and make sure all IBAMR tests pass (we do not run
   SAMRAI's old test suite and instead check for regressions in IBAMR)
2. Update `config/SAMRAI_config.h.in` to have the most recent version number.
   Nowadays we use the year, month, and day for the major, minor, and patchlevel
   numbers. To avoid problems with octal do not prefix numbers with `0` (e.g.,
   don't write `03` for the minor version number).
3. Create a tag with the same numbering scheme (e.g., 2026.10.21) as part of the
   most recent release. We tag commits on `master`.
