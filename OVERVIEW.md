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
  * `debug/`
    - `DFEAssert` - A method to throw an exception and halt a simulated DFE
  * `io/`
    - `AspectChangeIO` - save manager FIFOs by performing aspect change inside Kernels
    - `ShiftRegister` - perform serial -> parallel and parallel -> serial transformations via shift register
    - `VariableWidthIO` - a kernel input providing a variable number of words per cycle
    - `ZeroLatencyInput` - hide latency of input, allowing for data dependent control, e.g. run-length decoding
  * `mem/`
    - `BoxBuffer` - Buffer an N-dimensional stream of data in FMem and then read contiguous blocks out of it
    - `BoxBuffer1D` - Buffer a stream of data in FMem and read contiguous sections out of it
    - `ZeroLatencyMemory` - hide latency of FMem, allowing values to be read back on the next cycle
  * `pipeline/`
    - `FanoutLimiter` - Creates a tree of registers to minimise the fan-out of a KernelObect and improve timing
  * `sort/`
    - `BitonicSort`  - Sorts a list of DFEVars using a bitonic sort
  - `LargeStreamOffset` - large, negative stream offset backed by LMem to save FMem
* `lmem/`
  * `cpuaccess/` - simplified access to LMem from CPU
  - `MultiDimensionalAddressGenerator` - memory address generator in N dimensions
* `manager/`
  * `superfifo/` - extremely deep FIFO utilising LMem when necessary
  - `RoundRobin` - round robin data between outputs
* `network/`
  * `tcp/` - TCP framing (turn continuous TCP data into discrete frames)
* `statemachine/`
  * `collections/` - list, queue and stack implementations for state machines
* `utils/`
  - `TreeReduce` - reduce latency and save resources by reducing via a tree

