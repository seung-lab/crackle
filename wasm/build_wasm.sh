#!/bin/bash -xve

compile_options=(
    -Ithird_party/crc
    crackle_wasm.cc
     -O3
     -DNDEBUG
     --no-entry
     # -fno-exceptions
     -fno-rtti
     -s FILESYSTEM=0
     -s ALLOW_MEMORY_GROWTH=1 
     -s TOTAL_STACK=32768
     -s TOTAL_MEMORY=64kb
     -s EXPORTED_FUNCTIONS='["_crackle_compress","_crackle_decompress","_malloc","_free"]'
     -s MALLOC=emmalloc
     -s ENVIRONMENT=worker
     -s STANDALONE_WASM=1
     # -s ASSERTIONS=1 
     # -s SAFE_HEAP=1
     -std=c++2a
     -o libcrackle.wasm
)
em++ ${compile_options[@]}
