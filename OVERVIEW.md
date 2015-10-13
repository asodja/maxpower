Overview
========

A (non-exhaustive) list of the contents of the MaxPower library.

`maxpower/`
* `blas/l3/` - building blocks for dense matrix multiplication
* `hash/` - hash maps backed in FMem, LMem and QMem
* `kernel/`
  * `arithmetic/`
    - `ConstDenominator` - efficient division and modulus division by small, constant denominators
    - `FloatingPointAccumulator` - several floating point accumulator implementations
    - `FloatingPointMultiAdder` - save resources when adding three or more floating point numbers
  * `io/`
    - `AspectChangeIO` - save manager FIFOs by performing aspect change inside Kernels
    - `ZeroLatencyInput` - hide latency of input, allowing for data dependent control, e.g. run-length decoding
  * `mem/`
    - `BoxBuffer` - Buffer an N-dimensional stream of data in FMem and then read contiguous blocks out of it
    - `BoxBuffer1D` - Buffer a stream of data in FMem and read contiguous sections out of it
    - `ZeroLatencyMemory` - hide latency of FMem, allowing values to be read back on the next cycle
  * `LargeStreamOffset` - large, negative stream offset backed by LMem to save FMem
  * `TreeReduce` - reduce latency and save resources by reducing via a binary tree
* `lmem/`
  * `addressgenerators/` - memory address generators
  * `cpu_access/` - simplified access to LMem from CPU
  * `superfifo/` - extremely deep FIFO utilising LMem when necessary
* `manager/`
  - `RoundRobin` - round robin data between outputs
* `network/`
  * `tcp/` - TCP framing (turn continuous TCP data into discrete frames)
* `statemachine/`
  * `collections/` - list, queue and stack implementations for state machines
