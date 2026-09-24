# TH

CHICKEN Scheme bindings for the TH (Torch) tensor library
(`THTensor`, `THGenerator`, `THLongStorage`).

## Contents

* `th.scm` — the `th` module (compiled extension).
* `th.egg` — egg metadata / build options (includes a custom build step).
* `build-th` — custom build wrapper used by `chicken-install`.
* `tests/run.scm` — test suite.
* `TH/` — the vendored TH library sources (built automatically by the egg).

## Building and installing this egg

TH is bundled, so there is no separate TH installation step:
`chicken-install` compiles the vendored TH sources into a static archive
(`make -C TH static`) and links it directly into the extension. You only
need a C toolchain (`make`, `cc`, `ar`) and a BLAS library that provides
`libopenblas` (resolved via `ldconfig` at build time). Nothing is left to
find or install at runtime; the extension is self-contained.

If `/usr/local` is writable, `chicken-install` in this directory is
enough. Otherwise redirect the installation:

```
CHICKEN_INSTALL_PREFIX=$HOME/.local/lib/chicken \
CHICKEN_INSTALL_REPOSITORY=$HOME/.local/lib/chicken/12 \
chicken-install
```

CHICKEN 6 keeps its default repository under `/usr/local/lib/chicken/12`
(the old `…/11` was CHICKEN 5). If you installed CHICKEN elsewhere,
`chicken-install` reports the right value via the `###STATUS` above, or
check with `csi -p '(repository-path)'`.

Optional build tweaks are passed through `make` variables inside
`build-th` (e.g. `BLAS=`, `OMP=`, `CC=`) if defaults don't fit your
machine.

Note: `chicken-install -test` cannot validate this egg when the
repository is redirected (its test runner looks in the default
repository). Run the tests manually instead (below).

## Using it

After install, put the repository on the module path and import:

```
CHICKEN_REPOSITORY=$HOME/.local/lib/chicken/12 csi
# (import th)
```

Quick example:

```scheme
(import th)

(define a (list->float-tensor '(2 2) '(1.0 2.0 3.0 4.0)))
(define b (list->float-tensor '(2 2) '(5.0 6.0 7.0 8.0)))
(float-tensor->list (float-tensor-mm a b))    ; => (19.0 22.0 43.0 50.0)

(define t (make-float-tensor '(2 3)))
(float-tensor-fill! t 1.0)
(float-tensor-set2d! t 1 2 7.5)
(float-tensor-get2d t 1 2)                    ; => 7.5
(float-tensor-sum t)                          ; => 13.5

(define g (make-generator))
(generator-manual-seed! g 42)
(define r (make-float-tensor '(1000)))
(float-tensor-uniform! r g 0.0 1.0)           ; deterministic!
```

## API

There are seven tensor types; every operation is replicated per type
with the prefix `<type>-tensor-` where `<type>` is `float`, `double`,
`byte`, `char`, `short`, `int`, or `long`.

* Construction: `make-float-tensor`…`make-long-tensor`,
  `<type>-tensor-clone` (ranks 1–4; `make-*` raises for other ranks).
* Shape / metadata: `dims`, `strides`, `dim`, `size`, `stride`,
  `numel`, `contiguous?`, `data` (raw TH pointer), plus `tensor?`,
  `tensor->string` (numpy-format text), `tensor->pointer`.
* Element access: `ref`, `set!`, `get1d..get4d`, `set1d!..set4d!`
  (all validated; wrong rank/bounds raise a Scheme error, see Errors).
* Conversions: `<type>-tensor->list`, `list-><type>-tensor` (takes
  `(dims ls)`).
* In-place elementwise: `fill!`, `zero!`, `add!`, `mul!`, `div!`,
  `pow!`, `neg!`, `cinv!`, `abs!`, `sign!`, `sqrt!`, `cadd!`, `cmul!`,
  `cdiv!`, `cpow!`, `copy!` (`copy!` requires equal element counts).
* Reductions: `sum`, `mean`, `var`, `std`, `norm`, `min`, `max`,
  `median`, `dot`.
* Linear algebra: `addmm!`, `mm`.
* Generator: `make-generator`, `generator-valid?`, `generator-random`,
  `generator-seed!`, `generator-manual-seed!`, `generator-initial-seed`,
  `generator-uniform`, `generator-normal`.
* Random tensor fills: `random!`, `rand!`, `randn!`, `uniform!`,
  `normal!`, `exponential!`, `bernoulli!`, `geometric!`, `poisson!`,
  `cauchy!`.

### Integer type differences

TH only implements part of the math API for integer tensors, so these
names are missing for the integer types:

* `pow!`, `sqrt!` — float/double only.
* `mean`, `var`, `std`, `norm` — float/double only.
* `abs!` — short/int/long only (no `abs` for byte/char).
* `rand!`, `randn!`, `uniform!`, `normal!`, `exponential!`, `cauchy!`
  — float/double only. Integer random fills are `random!`,
  `bernoulli!`, `geometric!`, and `poisson!` (all require the generator
  as first argument, e.g. `(int-tensor-random! t g)`).

Element semantics follow C: `byte` stores the value modulo 256 and
roundtrips as 0–255, while `char` is the C `signed char` and
roundtrips as −128…127. Integer reductions return C `int`-precision
values; `sum`/`dot` come back as inexact doubles (TH's accumulator
type). `min`/`max`/`median` return exact integer values.

## Memory management

Every tensor and generator is allocated by TH and freed automatically
by a CHICKEN GC finalizer when no longer referenced — no explicit
`free` calls (this avoids double-free hazards). Tensors may outlive the
Scheme value that created them (e.g. `mm`, `clone`).

## Errors

TH's default error handler prints a message and calls `exit(-1)`.
This binding validates rank, bounds, and shapes on the Scheme side
first, so common mistakes raise Scheme exceptions instead of killing
the process. Untested/unsupported C-level misuse can still trip TH's
own checks and terminate the program.

## Running the tests

```
chicken-install            # build + install into your redirected repo
CHICKEN_REPOSITORY=<install-repo> CHICKEN_INCLUDE_PATH=<install-repo> \
  csi -s tests/run.scm
```

Point both variables at the directory `chicken-install` populated (the
value you used for `CHICKEN_INSTALL_REPOSITORY`). Note `csi` will ignore
the compiled extension in favour of a `th.scm` sitting in the current
directory, so run the script from outside the egg's source tree (or
delete generated `th.*.so`/`th.scm` copies) if imports misbehave.
