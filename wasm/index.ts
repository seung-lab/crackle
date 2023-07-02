/**
 * @license
 * Copyright 2023 William Silvermsith
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

const crackleWasmDataUrl = './libcrackle.wasm';

type TypedArray = Int8Array | Uint8Array | Int16Array | Uint16Array | Int32Array | Uint32Array | Uint8ClampedArray | Float32Array | Float64Array;

const libraryEnv = {
  emscripten_notify_memory_growth: function () {},
  proc_exit: (code: number) => {
    throw `proc exit: ${code}`;
  },
};

let wasmModule:any|null = null;

async function loadCrackleModule () {
  if (wasmModule !== null) {
    return wasmModule;
  }

  const response = await fetch(crackleWasmDataUrl);
  const wasmCode = await response.arrayBuffer();
  const m = await WebAssembly.instantiate(wasmCode, {
    env: libraryEnv,
    wasi_snapshot_preview1: libraryEnv,
  });
  (m.instance.exports._initialize as Function)();
  wasmModule = m;
  return m;
}

// not a full implementation of read header, just the parts we need
function readHeader(buffer: Uint8Array) 
  : {sx:number,sy:number,sz:number,dataWidth:number} 
{
  // check for header "crkl"
  const magic = (
       buffer[0] === 'c'.charCodeAt(0) && buffer[1] === 'r'.charCodeAt(0)
    && buffer[2] === 'k'.charCodeAt(0) && buffer[3] === 'l'.charCodeAt(0)
  );
  if (!magic) {
    throw new Error("crackle: didn't match magic numbers");
  }
  const format = buffer[4];
  if (format > 0) {
    throw new Error("crackle: didn't match format version");
  }

  const bufview = new DataView(buffer.buffer, 0);
  
  const format_bytes = bufview.getUint16(5, /*littleEndian=*/true);
  const dataWidth = Math.pow(2, format_bytes & 0b11);
  const sx = bufview.getUint32(7, /*littleEndian=*/true);
  const sy = bufview.getUint32(11, /*littleEndian=*/true);
  const sz = bufview.getUint32(15, /*littleEndian=*/true);

  return {sx,sy,sz,dataWidth};
}

function arrayType (dataWidth:number) : any {
  if (dataWidth === 1) {
    return Uint8Array;
  }
  else if (dataWidth === 2) {
    return Uint16Array;
  }
  else if (dataWidth === 4) {
    return Uint32Array;
  }
  else if (dataWidth === 8) {
    return BigUint64Array;
  }
}

export async function compressCrackle(
  buffer: Uint8Array,
  dataWidth:number,
  sx:number, sy:number, sz:number
) : Promise<Uint8Array> {
  const m = await loadCrackleModule();

  if (buffer.byteLength === 0) {
    throw new Error("crackle: Empty data buffer.");
  }
  else if (sx * sy * sz === 0) {
    throw new Error(`crackle: Invalid dimensions: <${sx},${sy},${sz}>`);
  }

  // heap must be referenced after creating bufPtr and imagePtr because
  // memory growth can detatch the buffer.
  let bufPtr = (m.instance.exports.malloc as Function)(buffer.byteLength);  
  let streamPtr = (m.instance.exports.malloc as Function)(buffer.byteLength);
  let heap = new Uint8Array((m.instance.exports.memory as WebAssembly.Memory).buffer);
  heap.set(buffer, bufPtr);

  const streamSize = (m.instance.exports.crackle_compress as Function)(
    bufPtr, dataWidth, 
    sx, sy, sz,
    streamPtr, buffer.byteLength
  );

  try {
    if (streamSize <= 0) {
      throw new Error(`crackle: Failed to encode image. encoder code: ${streamSize}`);
    }

    // Likewise, we reference memory.buffer instead of heap.buffer
    // because memory growth during decompress could have detached
    // the buffer.
    const stream = new Uint8Array(
      (m.instance.exports.memory as WebAssembly.Memory).buffer,
      streamPtr, streamSize
    );
    // copy the array so it can be memory managed by JS
    // and we can free the emscripten buffer
    return stream.slice(0);
  }
  finally {
    (m.instance.exports.free as Function)(bufPtr);
    (m.instance.exports.free as Function)(streamPtr);      
  }
}

export async function decompressCrackle(
  buffer: Uint8Array
) : Promise<Uint8Array> {

  const m = await loadCrackleModule();
  let {sx,sy,sz,dataWidth} = readHeader(buffer);

  const voxels = sx * sy * sz;
  const nbytes = voxels * dataWidth;
  if (nbytes < 0) {
    throw new Error(`crackle: Failed to decode image size. image size: ${nbytes}`);
  }

  // heap must be referenced after creating bufPtr and imagePtr because
  // memory growth can detatch the buffer.
  let bufPtr = (m.instance.exports.malloc as Function)(buffer.byteLength);
  const imagePtr = (m.instance.exports.malloc as Function)(nbytes);
  let heap = new Uint8Array((m.instance.exports.memory as WebAssembly.Memory).buffer);
  heap.set(buffer, bufPtr);

  const code = (m.instance.exports.crackle_decompress as Function)(
    bufPtr, buffer.byteLength, imagePtr, nbytes
  );

  try {
    if (code !== 0) {
      throw new Error(`crackle: Failed to decode image. decoder code: ${code}`);
    }

    // Likewise, we reference memory.buffer instead of heap.buffer
    // because memory growth during decompress could have detached
    // the buffer.
    let image = new Uint8Array(
      (m.instance.exports.memory as WebAssembly.Memory).buffer,
      imagePtr, nbytes
    );
    // copy the array so it can be memory managed by JS
    // and we can free the emscripten buffer
    image = image.slice(0);

    let ArrayType = arrayType(dataWidth);
    return new ArrayType(image.buffer);
  }
  finally {
    (m.instance.exports.free as Function)(bufPtr);
    (m.instance.exports.free as Function)(imagePtr);      
  }
}