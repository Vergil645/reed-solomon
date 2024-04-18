# Reed-Solomon

Repository contains efficient C implementation of Reed-Solomon erasure correction codes over $GF\left(2^{16}\right)$.

Encoder and decoder works with *symbol* sequences, where _symbol_ is a byte array of **even length**.

Encoder generates redundancy symbols for a given data sequence without changing it.
If you generate `r` redundancy symbols for `k` source symbols, you will be able to recover any `r` symbols from all `k + r` original symbols (source and redundancy).

Decoder recover missing symbols in a given sequence.

Main limitation is that `k + r` must be **less than 65536**.

## Performance compared to Random Linear Codes (RLC)

In the bellow you can find comparison of encoding and decoding time of implemented Reed-Solomon codes and RLC.

You can see that our implementation is **>2 times faster** than RLC.

| Compiler   | Optimizations level | CPU                             | Increase in coding speed relative to RLC | Increase in decoding speed relative to RLC |
|------------|:-------------------:|---------------------------------|:----------------------------------------:|:------------------------------------------:|
| gcc 12.2.0 | -O2                 | 12th Gen Intel® Core™ i7-12700F | 1.8                                      | 1.8                                        |
| gcc 12.2.0 | -O2                 | AMD Ryzen 5 3500U               | 1.7                                      | 1.6                                        |
| gcc 12.2.0 | -O3                 | 12th Gen Intel® Core™ i7-12700F | 2.6                                      | 2.5                                        |
| gcc 12.2.0 | -O3                 | AMD Ryzen 5 3500U               | 2.3                                      | 2.2                                        |

Also, our implementation consumes **60%** less energy than RLC. It was measured on 12th Gen Intel® Core™ i7-12700F CPU using `turbostat` utility.

## Interface and documentation

You can find documentation for encoder and decoder in [Reed-Solomon header](./include/rs/reed_solomon.h) file.

You can find documentation about symbol and sequence management in [memory include](./include/memory/) directory.

Also, you can generate HTML documentation using `Doxygen`:

```sh
doxygen
```

## Building library

### Debug

```sh
cmake -S. -Bbuild/debug -DBUILD_TESTING=true -DCMAKE_BUILD_TYPE=Debug
make -C build/debug
```

Running tests using `ctest`:

```sh
ctest --test-dir build/debug --stop-on-failure
```

Running `valgring` for test binaries:

```sh
find test/bin -name "test_*" | xargs -t -l1 valgrind -s -q --leak-check=full
```

### Release

```sh
cmake -S. -Bbuild/release -DBUILD_TESTING=false -DCMAKE_BUILD_TYPE=Release
make -C build/release
```

You can use option `-DADDITIONAL_C_FLAGS_RELEASE` to pass addition compiler options.

Example: `-DADDITIONAL_C_FLAGS_RELEASE="-O3"`.

### Output

Directory `lib` contains compiled library file (on Linux it is `librs.a`). You can link it with your project.

## Usage example

You can find example of using `librs` library in [example.c](./src/example.c) file.

If you want to include headers using `<...>`, like this

```c
#include <rs/reed_solomon.h>
```

you need to add `include` directory to include path of your compiler.
