/**
 * @license
 * Copyright 2021 William Silvermsith
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

#include "src/crackle.hpp"
#include <vector>

extern "C" {

int crackle_compress(
    unsigned char* buf, 
    unsigned int byte_width,
    unsigned int sx, unsigned int sy, unsigned int sz,
    unsigned char* out, unsigned int output_num_bytes
) {

    std::vector<unsigned char> outv;
    if (byte_width == 1) {
        outv = crackle::compress<uint8_t>(
            reinterpret_cast<uint8_t*>(buf), sx, sy, sz
        );
    }
    else if (byte_width == 2) {
        outv = crackle::compress<uint16_t>(
            reinterpret_cast<uint16_t*>(buf), sx, sy, sz
        );
    }
    else if (byte_width == 4) {
        outv = crackle::compress<uint32_t>(
            reinterpret_cast<uint32_t*>(buf), sx, sy, sz
        );
    }
    else if (byte_width == 8) {
        outv = crackle::compress<uint64_t>(
            reinterpret_cast<uint64_t*>(buf), sx, sy, sz
        );
    }

    unsigned int limit = std::min(
        static_cast<unsigned int>(outv.size()), output_num_bytes
    );

    for (unsigned int i = 0; i < limit; i++) {
        out[i] = outv[i];
    }

    return limit;
}

int crackle_decompress(
	unsigned char* buf, unsigned int num_bytes, 
    unsigned char* out, unsigned int output_num_bytes
) {
    crackle::decompress(buf, num_bytes, out, output_num_bytes);
    return 0;
}

}