/*
=== Kraken Decompressor for Windows ===
Copyright (C) 2016, Powzix

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "stdafx.h"
typedef unsigned char byte;
typedef unsigned char uint8;
typedef unsigned int uint32;
typedef unsigned __int64 uint64;
typedef signed __int64 int64;
typedef signed int int32;
typedef unsigned short uint16;

// Header in front of each 256k block
typedef struct KrakenHeader {
  // Type of decoder used, 6 means kraken
  int decoder_type;

  // Whether this block is uncompressed
  bool uncompressed;

  // Whether this block uses checksums.
  bool use_checksums;
} KrakenHeader;

// Additional header in front of each 256k block ("quantum").
typedef struct KrakenQuantumHeader {
  // The compressed size of this quantum. If this value is 0 it means
  // the quantum is a special quantum such as memset.
  uint32 compressed_size;
  // If checksums are enabled, holds the checksum.
  uint32 checksum;
  // Two flags
  uint8 flag1;
  uint8 flag2;
} KrakenQuantumHeader;

// Kraken decompression happens in two phases, first one decodes
// all the literals and copy lengths using huffman and second
// phase runs the copy loop. This holds the tables needed by stage 2.
typedef struct KrakenLzTable {
  // Stream of (literal, match) pairs. The flag byte contains
  // the length of the match, the length of the literal and whether
  // to use a recent distance.
  byte *flag_stream;
  int flag_stream_size;

  // Holds the actual distances in case we're not using a recent
  // distance.
  int *dist_stream;
  int dist_stream_size;

  // Holds the sequence of literals. All literal copying happens from
  // here.
  byte *lit_stream;
  int lit_stream_size;

  // Holds the lengths that do not fit in the flag stream. Both literal
  // lengths and match length are stored in the same array.
  int *len_stream;
  int len_stream_size;
} KrakenLzTable;


// Mermaid/Selkie decompression also happens in two phases, just like in Kraken,
// but the match copier works differently.
// Both Mermaid and Selkie use the same on-disk format, only the compressor
// differs.
typedef struct MermaidLzTable {
  // Flag stream. Format of flags:
  // Read flagbyte from |flag_stream|
  // If flagbyte >= 24:
  //   flagbyte & 0x80 == 0 : Read from |off16_stream| into |recent_distance|.
  //                   != 0 : Don't read offset.
  //   flagbyte & 7 = Number of literals to copy first from |lit_stream|.
  //   (flagbyte >> 3) & 0xF = Number of bytes to copy from |recent_distance|.
  //
  //  If flagbyte == 0 :
  //    Read byte L from |length_stream|
  //    If L > 251: L += 4 * Read word from |length_stream|
  //    L += 64
  //    Copy L bytes from |lit_stream|.
  //
  //  If flagbyte == 1 :
  //    Read byte L from |length_stream|
  //    If L > 251: L += 4 * Read word from |length_stream|
  //    L += 91
  //    Copy L bytes from match pointed by next offset from |off16_stream|
  //
  //  If flagbyte == 2 :
  //    Read byte L from |length_stream|
  //    If L > 251: L += 4 * Read word from |length_stream|
  //    L += 29
  //    Copy L bytes from match pointed by next offset from |off32_stream|, 
  //    relative to start of block.
  //    Then prefetch |off32_stream[3]|
  //
  //  If flagbyte > 2:
  //    L = flagbyte + 5
  //    Copy L bytes from match pointed by next offset from |off32_stream|,
  //    relative to start of block.
  //    Then prefetch |off32_stream[3]|
  const byte *flag_stream, *flag_stream_end;
  
  // Length stream
  const byte *length_stream;

  // Literal stream
  const byte *lit_stream, *lit_stream_end;

  // Near offsets
  const uint16 *off16_stream, *off16_stream_end;

  // Far offsets for current chunk
  uint32 *off32_stream, *off32_stream_end;
  
  // Holds the offsets for the two chunks
  uint32 *off32_stream_1, *off32_stream_2;
  uint32 off32_size_1, off32_size_2;

  // Flag offsets for next 64k chunk.
  uint32 flag_stream_2_offs, flag_stream_2_offs_end;
} MermaidLzTable;


typedef struct KrakenDecoder {
  // Updated after the |*_DecodeStep| function completes to hold
  // the number of bytes read and written.
  int src_used, dst_used;

  // Pointer to a 256k buffer that holds the intermediate state
  // in between decode phase 1 and 2.
  byte *phase_buf;
  int phase_buf_size;
} KrakenDecoder;

typedef struct HuffmanTree {
  // Number of symbols in the huffman tree. Always 8 for kraken.
  int num_syms;
  // Maximum allowed code length. Always 11 for kraken.
  int max_code_len;
  // Number of symbols actually used.
  int used_syms;
  // Index of last used symbol.
  int last_used_sym;
  // Min / Max code length
  int min_used_codelen, max_used_codelen;
  // If the huffman tree contains just one symbol this holds the symbol value,
  // then we can use memset and avoid storing anything.
  byte single_symbol;
  // Number of symbols for each code length
  int codelen_count[17];
  int F[17];
  int G[17];
  byte *symlen;
  byte *B;
  byte *C;
  // 1<<|max_code_length| entry mapping that maps a bit pattern to a code length.
  byte *bits_to_len;
  // 1<<|max_code_length| entry mapping that maps a bit pattern to a symbol.
  byte *bits_to_sym;
} HuffmanTree;

typedef struct BitReader {
  // |p| holds the current byte and |p_end| the end of the buffer.
  const byte *p, *p_end;
  // Bits accumulated so far
  uint32 bits;
  // Next byte will end up in the |bitpos| position in |bits|.
  int bitpos;
} BitReader;

typedef struct HuffReader {
  // Array to hold the output of the huffman read array operation
  byte *output, *output_end;
  // The code table
  HuffmanTree *huff;
  // We decode three parallel streams, two forwards, |src| and |src_mid|
  // while |src_end| is decoded backwards. 
  const byte *src, *src_mid, *src_end, *src_mid_org;
  int src_bitpos, src_mid_bitpos, src_end_bitpos;
  uint32 src_bits, src_mid_bits, src_end_bits;
} HuffReader;


// Allocate memory with a specific alignment
void *MallocAligned(size_t size, size_t alignment) {
  void *x = malloc(size + (alignment - 1) + sizeof(void*)), *x_org = x;
  if (x) {
    x = (void*)(((intptr_t)x + alignment - 1 + sizeof(void*)) & ~(alignment - 1));
    ((void**)x)[-1] = x_org;
  }
  return x;
}

// Free memory allocated through |MallocAligned|
void FreeAligned(void *p) {
  free(((void**)p)[-1]);
}

// Read more bytes to make sure we always have at least 24 bits in |bits|.
void BitReader_Refill(BitReader *bits) {
  while (bits->bitpos > 0) {
    bits->bits |= (bits->p < bits->p_end ? *bits->p : 0) << bits->bitpos;
    bits->bitpos -= 8;
    bits->p++;
  }
}

// Read more bytes to make sure we always have at least 24 bits in |bits|,
// used when reading backwards.
void BitReader_RefillBackwards(BitReader *bits) {
  while (bits->bitpos > 0) {
    bits->p--;
    bits->bits |= (bits->p >= bits->p_end ? *bits->p : 0) << bits->bitpos;
    bits->bitpos -= 8;
  }
}

// Refill bits then read a single bit.
int BitReader_ReadBit(BitReader *bits) {
  int r;
  BitReader_Refill(bits);
  r = bits->bits >> 31;
  bits->bits <<= 1;
  bits->bitpos += 1;
  return r;
}

// Read |n| bits without refilling.
int BitReader_ReadBitsNoRefill(BitReader *bits, int n) {
  int r = (bits->bits >> (32 - n));
  bits->bits <<= n;
  bits->bitpos += n;
  return r;
}

// Reads a gamma value.
// Assumes bitreader is already filled with at least 23 bits
int BitReader_ReadGamma(BitReader *bits) {
  unsigned long bitresult;
  int n;
  int r;
  if (bits->bits != 0) {
    _BitScanReverse(&bitresult, bits->bits);
    n = 31 - bitresult;
  } else {
    n = 32;
  }
  n = 2 * n + 2;
  assert(n < 24);
  bits->bitpos += n;
  r = bits->bits >> (32 - n);
  bits->bits <<= n;
  return r - 2;
}

// Reads a gamma value with |forced| number of forced bits.
int BitReader_ReadGammaX(BitReader *bits, int forced) {
  unsigned long bitresult;
  int n, r;
  if (bits->bits != 0) {
    _BitScanReverse(&bitresult, bits->bits);
    n = 31 - bitresult;
    assert(n < 24);
    r = (bits->bits >> (31 - n - forced)) + ((n - 1) << forced);
    bits->bits <<= n + forced + 1;
    bits->bitpos += n + forced + 1;
    return r;
  }
  return 0;
}

// Reads a distance code parametrized by |v|.
uint32 BitReader_ReadDistance(BitReader *bits, uint32 v) {
  uint32 w, m, n, rv;
  if (v < 0xF0) {
    n = (v >> 4) + 4;
    w = _rotl(bits->bits | 1, n);
    bits->bitpos += n;
    m = (2 << n) - 1;
    bits->bits = w & ~m;
    rv = ((w & m) << 4) + (v & 0xF) - 248;
  } else {
    n = v - 0xF0 + 4;
    w = _rotl(bits->bits | 1, n);
    bits->bitpos += n;
    m = (2 << n) - 1;
    bits->bits = w & ~m;
    rv = 8322816 + ((w & m) << 12);
    BitReader_Refill(bits);
    rv += (bits->bits >> 20);
    bits->bitpos += 12;
    bits->bits <<= 12;
  }
  BitReader_Refill(bits);
  return rv;
}

// Reads a distance code parametrized by |v|, backwards.
uint32 BitReader_ReadDistanceB(BitReader *bits, uint32 v) {
  uint32 w, m, n, rv;
  if (v < 0xF0) {
    n = (v >> 4) + 4;
    w = _rotl(bits->bits | 1, n);
    bits->bitpos += n;
    m = (2 << n) - 1;
    bits->bits = w & ~m;
    rv = ((w & m) << 4) + (v & 0xF) - 248;
  } else {
    n = v - 0xF0 + 4;
    w = _rotl(bits->bits | 1, n);
    bits->bitpos += n;
    m = (2 << n) - 1;
    bits->bits = w & ~m;
    rv = 8322816 + ((w & m) << 12);
    BitReader_RefillBackwards(bits);
    rv += (bits->bits >> (32 - 12));
    bits->bitpos += 12;
    bits->bits <<= 12;
  }
  BitReader_RefillBackwards(bits);
  return rv;
}

// Reads a length code.
bool BitReader_ReadLength(BitReader *bits, uint32 *v) {
  unsigned long bitresult;
  int n;
  uint32 rv;
  _BitScanReverse(&bitresult, bits->bits);
  n = 31 - bitresult;
  if (n > 12) return false;
  bits->bitpos += n;
  bits->bits <<= n;
  BitReader_Refill(bits);
  n += 7;
  bits->bitpos += n;
  rv = (bits->bits >> (32 - n)) - 64;
  bits->bits <<= n;
  *v = rv;
  BitReader_Refill(bits);
  return true;
}

// Reads a length code, backwards.
bool BitReader_ReadLengthB(BitReader *bits, uint32 *v) {
  unsigned long bitresult;
  int n;
  uint32 rv;
  _BitScanReverse(&bitresult, bits->bits);
  n = 31 - bitresult;
  if (n > 12) return false;
  bits->bitpos += n;
  bits->bits <<= n;
  BitReader_RefillBackwards(bits);
  n += 7;
  bits->bitpos += n;
  rv = (bits->bits >> (32 - n)) - 64;
  bits->bits <<= n;
  *v = rv;
  BitReader_RefillBackwards(bits);
  return true;
}

int Log2RoundUp(uint32 v) {
  if (v > 1) {
    unsigned long idx;
    _BitScanReverse(&idx, v - 1);
    return idx + 1;
  } else {
    return 0;
  }
}

#define ALIGN_16(x) (((x)+15)&~15)
#define COPY_64(d, s) {*(uint64*)(d) = *(uint64*)(s); }
#define COPY_64_ADD(d, s, t) _mm_storel_epi64((__m128i *)(d), _mm_add_epi8(_mm_loadl_epi64((__m128i *)(s)), _mm_loadl_epi64((__m128i *)(t))))


// Allocate a huffman tree.
HuffmanTree *HuffmanTree_Allocate(int num_syms, int max_code_len, void *memory) {
  HuffmanTree *huff = (HuffmanTree*)memory;
  memset(huff, 0, sizeof(HuffmanTree));
  huff->num_syms = num_syms;
  huff->max_code_len = max_code_len;
  huff->symlen = (byte*)(huff + 1);
  huff->B = huff->symlen + ALIGN_16(num_syms + 1);
  huff->C = (byte*)huff->B;
  huff->bits_to_len = (byte*)huff->C + ALIGN_16(2 * (num_syms + 1));
  huff->bits_to_sym = (byte*)(huff->bits_to_len + (1 << max_code_len) + 16);
  return huff;
}

// Read and decode the huffman tree from the bit reader |bits|.
bool HuffmanTree_Decode(HuffmanTree *huff, BitReader *bits) {
  if (BitReader_ReadBit(bits)) {
    int num_syms = huff->num_syms;
    int symbits = Log2RoundUp(num_syms);
    int n, sym = 0, codelen;
    int avg_bits_x4 = symbits * 4;
    int min_codelen = 32, max_codelen = 0;
    unsigned int p;
    int forced_bits = BitReader_ReadBitsNoRefill(bits, 2);
    if (BitReader_ReadBit(bits))
      goto SKIP_INITIAL_ZEROS;
    for (;;) {
      // Run of zeros
      n = BitReader_ReadGamma(bits) + 1;
      if (sym + n > huff->num_syms)
        return false;
      memset(&huff->symlen[sym], 0, n);
      sym += n;
      if (sym >= huff->num_syms)
        break;
      BitReader_Refill(bits);
SKIP_INITIAL_ZEROS:
      // Then a run of non-zero
      n = BitReader_ReadGamma(bits) + 1;
      if (sym + n > huff->num_syms)
        return false;
      BitReader_Refill(bits);
      do {
        p = BitReader_ReadGammaX(bits, forced_bits);
        BitReader_Refill(bits);
        codelen = (-(int)(p & 1) ^ (p >> 1)) + ((avg_bits_x4 + 2) >> 2);
        if (codelen < 0 || codelen > 16)
          return false;
        huff->used_syms++;
        huff->codelen_count[codelen]++;

        if (codelen < min_codelen)
          min_codelen = codelen;
        
        if (codelen > max_codelen)
          max_codelen = codelen;

        avg_bits_x4 = codelen + ((3 * avg_bits_x4 + 2) >> 2);
        huff->last_used_sym = sym;
        huff->symlen[sym++] = codelen;
      } while (--n);
      if (sym >= huff->num_syms)
        break;
    }
    huff->min_used_codelen = min_codelen;
    huff->max_used_codelen = max_codelen;
    return true;
  } else {
    int sym_bits = Log2RoundUp(huff->num_syms);
    huff->used_syms = BitReader_ReadBitsNoRefill(bits, sym_bits);

    if (huff->used_syms > huff->num_syms)
      return false;

    if (huff->used_syms == 0)
      return true;

    if (huff->used_syms != 1) {
      int codelen_bits, sym, codelen;
      int i;
      int min_codelen = 32, max_codelen = 0;

      memset(huff->symlen, 0, huff->num_syms);

      codelen_bits = BitReader_ReadBitsNoRefill(bits, 3);
      if (codelen_bits > 5)
        return false;

      for (i = 0; i < huff->used_syms; i++) {
        BitReader_Refill(bits);

        sym = BitReader_ReadBitsNoRefill(bits, sym_bits);
        codelen = BitReader_ReadBitsNoRefill(bits, codelen_bits) + 1;

        if (codelen > 16 || sym >= huff->num_syms)
          return false;

        huff->codelen_count[codelen]++;
        huff->last_used_sym = sym;
        huff->symlen[sym] = codelen;

        if (codelen < min_codelen)
          min_codelen = codelen;

        if (codelen > max_codelen)
          max_codelen = codelen;
      }

      huff->min_used_codelen = min_codelen;
      huff->max_used_codelen = max_codelen;
      return true;
    } else {
      int s;
      BitReader_Refill(bits);
      huff->single_symbol = s = BitReader_ReadBitsNoRefill(bits, sym_bits);
      if (s >= huff->num_syms)
        return false;
      huff->last_used_sym = s;
      huff->symlen[s] = 0;
      huff->C[0] = s;
      return true;
    }
  }
}

void HuffmanTree_AssignSyms(byte *symlen, int num_used_syms, int *Z, uint8 *C) {
  int i;
  for (i = 0; i < num_used_syms; i++)
    C[Z[symlen[i]]++] = i;
}

bool HuffmanTree_SetupBits(HuffmanTree *huff) {
  int i, i_end;
  int local_arr[17];
  int x, y, z;

  assert(huff->used_syms > 1);

  if (huff->min_used_codelen == 0 || huff->max_used_codelen == 0 || huff->max_used_codelen > huff->max_code_len)
    return false;

  memset(huff->F, 0, sizeof(huff->F[0]) * huff->min_used_codelen);

  x = y = z = 0;
  for (i = huff->min_used_codelen, i_end = huff->max_used_codelen; i <= i_end; i++) {
    local_arr[i] = x;
    z = huff->codelen_count[i] + y;
    huff->G[i] = y - x;
    huff->F[i] = z << (32 - i);
    y = 2 * z;
    x += huff->codelen_count[i];
  }
  if (z != (1 << huff->max_used_codelen))
    return false;

  memset(&huff->F[huff->max_used_codelen], 0xff, sizeof(huff->F[0]) * (32 - huff->max_used_codelen));
  local_arr[0] = huff->used_syms;
  HuffmanTree_AssignSyms(huff->symlen, huff->last_used_sym + 1, local_arr, huff->C);
  return true;
}

bool HuffmanTree_MakeLookupTables(HuffmanTree *huff) {
  int i, shift, m, n, p;
  uint8 *cursym;
  if (!HuffmanTree_SetupBits(huff))
    return false;
  
  assert(huff->used_syms > 1);

  for (i = huff->min_used_codelen, cursym = huff->C, p = 0; i != huff->max_code_len; i++) {
    n = huff->codelen_count[i];
    if (n) {
      shift = huff->max_code_len - i;
      memset(huff->bits_to_len + p, i, n << shift);
      m = 1 << shift;
      do {
        memset(huff->bits_to_sym + p, *cursym++, m);
        p += m;
      } while (--n);
    }
  }
  if ((n = huff->codelen_count[i]) != 0) {
    memset(huff->bits_to_len + p, i, n);
    memmove(huff->bits_to_sym + p, cursym, sizeof(cursym[0]) * n);
  }
  return true;
}

KrakenDecoder *Kraken_Create() {
  int phase_buf_size = 0x40000;
  int memory_needed = sizeof(KrakenDecoder) + phase_buf_size;
  KrakenDecoder *dec = (KrakenDecoder*)MallocAligned(memory_needed, 16);
  memset(dec, 0, sizeof(KrakenDecoder));
  dec->phase_buf_size = phase_buf_size;
  dec->phase_buf = (byte*)(dec + 1);
  return dec;
}

void Kraken_Destroy(KrakenDecoder *kraken) {
  FreeAligned(kraken);
}

const byte *Kraken_ParseHeader(KrakenHeader *hdr, const byte *p) {
  int b = p[0];
  if ((b & 0xF) != 0xC) return NULL;
  if (((b >> 4) & 3) != 0) return NULL;
  hdr->uncompressed = (b >> 6) & 1;
  b = p[1];
  hdr->decoder_type = b & 0x3F;
  hdr->use_checksums = !!(b >> 7);
  if (hdr->decoder_type != 6 && hdr->decoder_type != 10) return NULL;
  return p + 2;
}

const byte *Kraken_ParseQuantumHeader(KrakenQuantumHeader *hdr, const byte *p, bool use_checksum) {
  uint32 v = (p[0] << 16) | (p[1] << 8) | p[2];
  uint32 size = v & 0x3FFFF;
  if (size != 0x3ffff) {
    hdr->compressed_size = size + 1;
    hdr->flag1 = (v >> 18) & 1;
    hdr->flag2 = (v >> 19) & 1;
    if (use_checksum) {
      hdr->checksum = (p[3] << 16) | (p[4] << 8) | p[5];
      return p + 6;
    } else {
      return p + 3;
    }
  }
  v >>= 18;
  if (v == 1) {
    hdr->checksum = p[3];
    hdr->compressed_size = 0;
    return p + 4;
  }
  return NULL;

}

uint32 Kraken_GetCrc(const byte *p, size_t p_size) {
  // TODO: implement
  return 0;
}

// Rearranges elements in the input array so that bits in the index
// get flipped.
static void ReverseBitsArray2048(const byte *input, byte *output) {
  static const uint8 offsets[32] = {
    0,    0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
    0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8
  };
  __m128i t0, t1, t2, t3, s0, s1, s2, s3;
  int i, j;
  for(i = 0; i != 32; i++) {
    j = offsets[i];
    t0 = _mm_unpacklo_epi8(
      _mm_loadl_epi64((const __m128i *)&input[j]),
      _mm_loadl_epi64((const __m128i *)&input[j + 256]));
    t1 = _mm_unpacklo_epi8(
      _mm_loadl_epi64((const __m128i *)&input[j + 512]),
      _mm_loadl_epi64((const __m128i *)&input[j + 768]));
    t2 = _mm_unpacklo_epi8(
      _mm_loadl_epi64((const __m128i *)&input[j + 1024]),
      _mm_loadl_epi64((const __m128i *)&input[j + 1280]));
    t3 = _mm_unpacklo_epi8(
      _mm_loadl_epi64((const __m128i *)&input[j + 1536]),
      _mm_loadl_epi64((const __m128i *)&input[j + 1792]));

    s0 = _mm_unpacklo_epi8(t0, t1);
    s1 = _mm_unpacklo_epi8(t2, t3);
    s2 = _mm_unpackhi_epi8(t0, t1);
    s3 = _mm_unpackhi_epi8(t2, t3);

    t0 = _mm_unpacklo_epi8(s0, s1);
    t1 = _mm_unpacklo_epi8(s2, s3);
    t2 = _mm_unpackhi_epi8(s0, s1);
    t3 = _mm_unpackhi_epi8(s2, s3);

    _mm_storel_epi64((__m128i *)&output[0], t0);
    _mm_storeh_pi((__m64*)&output[1024], _mm_castsi128_ps(t0));
    _mm_storel_epi64((__m128i *)&output[256], t1);
    _mm_storeh_pi((__m64*)&output[1280], _mm_castsi128_ps(t1));
    _mm_storel_epi64((__m128i *)&output[512], t2);
    _mm_storeh_pi((__m64*)&output[1536], _mm_castsi128_ps(t2));
    _mm_storel_epi64((__m128i *)&output[768], t3);
    _mm_storeh_pi((__m64*)&output[1792], _mm_castsi128_ps(t3));
    output += 8;
  }
}

bool Kraken_DecodeBytesCore(HuffReader *hr) {
  // Create the reversed arrays
  uint8 bits2len_rev[2048];
  uint8 bits2sym_rev[2048];

  ReverseBitsArray2048(hr->huff->bits_to_len, bits2len_rev);
  ReverseBitsArray2048(hr->huff->bits_to_sym, bits2sym_rev);

  const byte *src = hr->src;
  uint32 src_bits = hr->src_bits;
  int src_bitpos = hr->src_bitpos;

  const byte *src_mid = hr->src_mid;
  uint32 src_mid_bits = hr->src_mid_bits;
  int src_mid_bitpos = hr->src_mid_bitpos;

  const byte *src_end = hr->src_end;
  uint32 src_end_bits = hr->src_end_bits;
  int src_end_bitpos = hr->src_end_bitpos;

  int k, n;

  byte *dst = hr->output;
  byte *dst_end = hr->output_end;

  if (src > src_mid)
    return false;

  if (hr->src_end - src_mid >= 4 && dst_end - dst >= 6) {
    dst_end -= 5;
    src_end -= 4;

    while (dst < dst_end && src <= src_mid && src_mid <= src_end) {
      src_bits |= *(uint32*)src << src_bitpos;
      src += (31 - src_bitpos) >> 3;

      src_end_bits |= _byteswap_ulong(*(uint32*)src_end) << src_end_bitpos;
      src_end -= (31 - src_end_bitpos) >> 3;

      src_mid_bits |= *(uint32*)src_mid << src_mid_bitpos;
      src_mid += (31 - src_mid_bitpos) >> 3;

      src_bitpos |= 0x18;
      src_end_bitpos |= 0x18;
      src_mid_bitpos |= 0x18;

      k = src_bits & 0x7FF;
      n = bits2len_rev[k];
      src_bits >>= n;
      src_bitpos -= n;
      dst[0] = bits2sym_rev[k];

      k = src_end_bits & 0x7FF;
      n = bits2len_rev[k];
      src_end_bits >>= n;
      src_end_bitpos -= n;
      dst[1] = bits2sym_rev[k];

      k = src_mid_bits & 0x7FF;
      n = bits2len_rev[k];
      src_mid_bits >>= n;
      src_mid_bitpos -= n;
      dst[2] = bits2sym_rev[k];

      k = src_bits & 0x7FF;
      n = bits2len_rev[k];
      src_bits >>= n;
      src_bitpos -= n;
      dst[3] = bits2sym_rev[k];

      k = src_end_bits & 0x7FF;
      n = bits2len_rev[k];
      src_end_bits >>= n;
      src_end_bitpos -= n;
      dst[4] = bits2sym_rev[k];

      k = src_mid_bits & 0x7FF;
      n = bits2len_rev[k];
      src_mid_bits >>= n;
      src_mid_bitpos -= n;
      dst[5] = bits2sym_rev[k];
      dst += 6;
    }
    dst_end += 5;

    src -= src_bitpos >> 3;
    src_bitpos &= 7;

    src_end += 4 + (src_end_bitpos >> 3);
    src_end_bitpos &= 7;

    src_mid -= src_mid_bitpos >> 3;
    src_mid_bitpos &= 7;
  }
  for(;;) {
    if (dst >= dst_end)
      break;

    if (src_mid - src <= 1) {
      if (src_mid - src == 1)
        src_bits |= *src << src_bitpos;
    } else {
      src_bits |= *(uint16 *)src << src_bitpos;
    }
    k = src_bits & 0x7FF;
    n = bits2len_rev[k];
    src_bitpos -= n;
    src_bits >>= n;
    *dst++ = bits2sym_rev[k];
    src += (7 - src_bitpos) >> 3;
    src_bitpos &= 7;

    if (dst < dst_end) {
      if (src_end - src_mid <= 1) {
        if (src_end - src_mid == 1) {
          src_end_bits |= *src_mid << src_end_bitpos;
          src_mid_bits |= *src_mid << src_mid_bitpos;
        }
      } else {
        unsigned int v = *(uint16*)(src_end - 2);
        src_end_bits |= (((v >> 8) | (v << 8)) & 0xffff) << src_end_bitpos;
        src_mid_bits |= *(uint16*)src_mid << src_mid_bitpos;
      }
      n = bits2len_rev[src_end_bits & 0x7FF];
      *dst++ = bits2sym_rev[src_end_bits & 0x7FF];
      src_end_bitpos -= n;
      src_end_bits >>= n;
      src_end -= (7 - src_end_bitpos) >> 3;
      src_end_bitpos &= 7;
      if (dst < dst_end) {
        n = bits2len_rev[src_mid_bits & 0x7FF];
        *dst++ = bits2sym_rev[src_mid_bits & 0x7FF];
        src_mid_bitpos -= n;
        src_mid_bits >>= n;
        src_mid += (7 - src_mid_bitpos) >> 3;
        src_mid_bitpos &= 7;
      }
    }
    if (src > src_mid || src_mid > src_end)
      return false;
  }
  if (src != hr->src_mid_org || src_end != src_mid)
    return false;
  return true;
}

bool Kraken_DecodeBytes_Type12(const byte *src, const byte *src_end, byte *output, int output_size, int type) {
  // Temporary storage for huffman table
  byte temp_mem[8192];
  HuffmanTree *huff;
  BitReader bits;
  int bytes_parsed, half_output_size;
  uint32 split_left, split_mid, split_right;
  const byte *src_mid;
  HuffReader hr;
  huff = HuffmanTree_Allocate(256, 11, temp_mem);
  bits.bitpos = 24;
  bits.bits = 0;
  bits.p = src;
  bits.p_end = src_end;
  if (BitReader_ReadBit(&bits) ||
      !HuffmanTree_Decode(huff, &bits))
    return false;
  bytes_parsed = bits.p - src - ((24 - bits.bitpos) / 8);
  src += bytes_parsed;
  if (huff->used_syms <= 1) {
    if (src != src_end)
      return false;
    memset(output, huff->single_symbol, output_size);
    return true;
  }
  if (!HuffmanTree_MakeLookupTables(huff))
    return false;

  if (type == 1) {
    if (src + 3 > src_end)
      return false;
    split_mid = *(uint16*)src;
    src += 2;
    hr.huff = huff;
    hr.output = output;
    hr.output_end = output + output_size;
    hr.src = src;
    hr.src_end = src_end;
    hr.src_mid_org = hr.src_mid = src + split_mid;
    hr.src_bitpos = 0;
    hr.src_bits = 0;
    hr.src_mid_bitpos = 0;
    hr.src_mid_bits = 0;
    hr.src_end_bitpos = 0;
    hr.src_end_bits = 0;
    return Kraken_DecodeBytesCore(&hr);
  } else {
    if (src + 6 > src_end)
      return false;

    half_output_size = (output_size + 1) >> 1;
    split_mid = *(uint32*)src & 0xFFFFFF;
    src += 3;
    if (split_mid > (src_end - src))
      return false;
    src_mid = src + split_mid;
    split_left = *(uint16*)src;
    src += 2;
    if (src_mid - src < split_left + 2 ||
        src_end - src_mid < 3)
      return false;
    split_right = *(uint16*)src_mid;
    if (src_end - (src_mid + 2) < split_right + 2)
      return false;

    hr.huff = huff;
    hr.output = output;
    hr.output_end = output + half_output_size;
    hr.src = src;
    hr.src_end = src_mid;
    hr.src_mid_org = hr.src_mid = src + split_left;
    hr.src_bitpos = 0;
    hr.src_bits = 0;
    hr.src_mid_bitpos = 0;
    hr.src_mid_bits = 0;
    hr.src_end_bitpos = 0;
    hr.src_end_bits = 0;
    if (!Kraken_DecodeBytesCore(&hr))
      return false;

    hr.huff = huff;
    hr.output = output + half_output_size;
    hr.output_end = output + output_size;
    hr.src = src_mid + 2;
    hr.src_end = src_end;
    hr.src_mid_org = hr.src_mid = src_mid + 2 + split_right;
    hr.src_bitpos = 0;
    hr.src_bits = 0;
    hr.src_mid_bitpos = 0;
    hr.src_mid_bits = 0;
    hr.src_end_bitpos = 0;
    hr.src_end_bits = 0;
    if (!Kraken_DecodeBytesCore(&hr))
      return false;

    return true;
  }
}

int Kraken_DecodeBytes(byte **output, const byte *src, const byte *src_end, int *decoded_size, int output_size) {
  uint32 lo_bits;
  int src_size, num_elems;
  int type;

  if (src_end - src < 5)
    return -1;
  type = src[0] >> 5;
  switch (type) {
  case 0:
    src_size = ((src[0] << 16) | (src[1] << 8) | src[2]) & 0x3FFFF;
    src += 3;
    if (src_size > output_size || src_end - src < src_size)
      return -1;
    *decoded_size = src_size;
    *output = (byte*)src;
    return 3 + src_size;
  case 1:
  case 2:
    lo_bits = _byteswap_ulong(*(int*)(src + 1));
    src += 5;
    src_size = lo_bits & 0x3FFFF;
    if (src_end - src < src_size)
      return -1;
    num_elems = (((lo_bits >> 18) | (src[-5] << 14)) & 0x3FFFF) + 1;
    if (num_elems > output_size || src_size >= num_elems)
      return -1;
    if (output != NULL && !Kraken_DecodeBytes_Type12(src, src + src_size, *output, num_elems, type))
      return -1;
    *decoded_size = num_elems;
    return 5 + src_size;
  }
  return -1;
}

// Unpacks the packed 8 bit offset and lengths into 32 bit.
bool Kraken_UnpackOffsets(const byte *src, const byte *src_end,
                          const byte *packed_dist_stream, int packed_dist_stream_size,
                          const byte *packed_litlen_stream, int packed_litlen_stream_size,
                          int *dist_stream, int *len_stream, byte *output_end) {


  BitReader bits_a, bits_b;
  unsigned long idx;
  int n, i;
  int extralong_litlen_stream_size;
  uint32 *extralong_litlen_stream;

  bits_a.bitpos = 24;
  bits_a.bits = 0;
  bits_a.p = src;
  bits_a.p_end = src_end;
  BitReader_Refill(&bits_a);

  bits_b.bitpos = 24;
  bits_b.bits = 0;
  bits_b.p = src_end;
  bits_b.p_end = src;
  BitReader_RefillBackwards(&bits_b);

  _BitScanReverse(&idx, bits_b.bits);
  n = 31 - idx;
  if (n > 18) return false;
  bits_b.bitpos += n;
  bits_b.bits <<= n;
  BitReader_RefillBackwards(&bits_b);
  n++;
  extralong_litlen_stream_size = (bits_b.bits >> (32 - n)) - 1;
  bits_b.bitpos += n;
  bits_b.bits <<= n;
  BitReader_RefillBackwards(&bits_b);

  for (i = 0; i < packed_dist_stream_size - 1; i += 2) {
    dist_stream[i + 0] = -(int32)BitReader_ReadDistance(&bits_a, packed_dist_stream[i + 0]);
    dist_stream[i + 1] = -(int32)BitReader_ReadDistanceB(&bits_b, packed_dist_stream[i + 1]);
  }
  if (i < packed_dist_stream_size)
    dist_stream[i + 0] = -(int32)BitReader_ReadDistance(&bits_a, packed_dist_stream[i + 0]);
  
  extralong_litlen_stream = (uint32*)((intptr_t)output_end & ~3) - extralong_litlen_stream_size;
  for (i = 0; i < extralong_litlen_stream_size - 1; i += 2) {
    if (!BitReader_ReadLength(&bits_a, &extralong_litlen_stream[i + 0])) return false;
    if (!BitReader_ReadLengthB(&bits_b, &extralong_litlen_stream[i + 1])) return false;
  }
  if (i < extralong_litlen_stream_size)
    BitReader_ReadLength(&bits_a, &extralong_litlen_stream[i + 0]);

  bits_a.p -= (24 - bits_a.bitpos) >> 3;
  bits_b.p += (24 - bits_b.bitpos) >> 3;

  if (bits_a.p != bits_b.p)
    return false;

  for (i = 0; i < packed_litlen_stream_size; i++) {
    uint32 v = packed_litlen_stream[i];
    if (v == 255)
      v = *extralong_litlen_stream++ + 255;
    len_stream[i] = v + 3;
  }
  return true;
}

bool Kraken_ReadLzTable(int mode,
                        const byte *src, const byte *src_end,
                        byte *dst, int dst_size, int offset,
                        byte *output, byte *output_end, KrakenLzTable *lztable) {
  byte *out;
  int decode_count, n;
  byte *packed_dist_stream, *packed_litlen_stream;

  if (mode > 1)
    return false;

  if (src_end - src < 13)
    return false;

  if (offset == 0) {
    COPY_64(dst, src);
    dst += 8;
    src += 8;
  }

  // Decode lit stream
  out = output;
  n = Kraken_DecodeBytes(&out, src, src_end, &decode_count, output_end - output);
  if (n < 0)
    return false;
  src += n;
  lztable->lit_stream = out;
  lztable->lit_stream_size = decode_count;
  output += decode_count;

  // Decode flag stream
  out = output;
  n = Kraken_DecodeBytes(&out, src, src_end, &decode_count, output_end - output);
  if (n < 0 || decode_count == 0)
    return false;
  src += n;
  lztable->flag_stream = out;
  lztable->flag_stream_size = decode_count;
  output += decode_count;

  // Decode packed distance stream
  packed_dist_stream = output;
  n = Kraken_DecodeBytes(&packed_dist_stream, src, src_end, &lztable->dist_stream_size, output_end - output);
  if (n < 0)
    return false;
  src += n;
  output += lztable->dist_stream_size;

  // Decode packed litlen stream
  packed_litlen_stream = output;
  n = Kraken_DecodeBytes(&packed_litlen_stream, src, src_end, &lztable->len_stream_size, output_end - output);
  if (n < 0)
    return false;
  src += n;
  output += lztable->len_stream_size;

  // Reserve memory for final dist stream
  output = (byte*)(((intptr_t)output + 3) & ~3);
  lztable->dist_stream = (int*)output;
  output += lztable->dist_stream_size * 4;

  // Reserve memory for final litlen stream
  output = (byte*)(((intptr_t)output + 15) & ~15);
  lztable->len_stream = (int*)output;
  output += lztable->len_stream_size * 4;
  if (output > output_end)
    return false;

  return Kraken_UnpackOffsets(src, src_end, packed_dist_stream, lztable->dist_stream_size,
                            packed_litlen_stream, lztable->len_stream_size,
                            lztable->dist_stream, lztable->len_stream, output_end);
}


// Note: may access memory out of bounds on invalid input.
bool Kraken_ProcessLzRuns_Type0(KrakenLzTable *lzt, byte *dst, byte *dst_end, byte *dst_start) {
  const byte *flag_stream = lzt->flag_stream,
    *flag_stream_end = flag_stream + lzt->flag_stream_size;
  const int *len_stream = lzt->len_stream;
  const int *litlen_stream_end = lzt->len_stream + lzt->len_stream_size;
  const byte *lit_stream = lzt->lit_stream;
  const byte *lit_stream_end = lzt->lit_stream + lzt->lit_stream_size;
  const int *distance_stream = lzt->dist_stream;
  const int *distance_stream_end = lzt->dist_stream + lzt->dist_stream_size;
  const byte *copyfrom;
  uint32 final_len;
  int32 distance;
  int32 recent_distance[7];
  int32 last_lit_dist;

  recent_distance[3] = -8;
  recent_distance[4] = -8;
  recent_distance[5] = -8;
  last_lit_dist = -8;

  while (flag_stream < flag_stream_end) {
    uint32 f = *flag_stream++;
    uint32 litlen = f & 3;
    uint32 disttype = f >> 6;
    uint32 copy_len = (f >> 2) & 0xF;

    // use cmov
    uint32 next_long_length = *len_stream;
    const int *next_len_stream = len_stream + 1;

    len_stream = (litlen == 3) ? next_len_stream : len_stream;
    litlen = (litlen == 3) ? next_long_length : litlen;
    recent_distance[6] = *distance_stream;

    COPY_64_ADD(dst, lit_stream, &dst[last_lit_dist]);
    if (litlen > 8) {
      COPY_64_ADD(dst + 8, lit_stream + 8, &dst[last_lit_dist + 8]);
      if (litlen > 16) {
        COPY_64_ADD(dst + 16, lit_stream + 16, &dst[last_lit_dist + 16]);
        if (litlen > 24) {
          do {
            COPY_64_ADD(dst + 24, lit_stream + 24, &dst[last_lit_dist + 24]);
            litlen -= 8;
            dst += 8;
            lit_stream += 8;
          } while (litlen > 24);
        }
      }
    }
    dst += litlen;
    lit_stream += litlen;

    distance = recent_distance[disttype + 3];
    recent_distance[disttype + 3] = recent_distance[disttype + 2];
    recent_distance[disttype + 2] = recent_distance[disttype + 1];
    recent_distance[disttype + 1] = recent_distance[disttype + 0];
    recent_distance[3] = distance;
    last_lit_dist = distance;

    distance_stream = (int*)((intptr_t)distance_stream + ((disttype + 1) & 4));

    if ((uintptr_t)distance < (uintptr_t)(dst_start - dst))
      return false; // distance out of bounds

    copyfrom = dst + distance;
    if (copy_len != 15) {
      COPY_64(dst, copyfrom);
      COPY_64(dst + 8, copyfrom + 8);
      dst += copy_len + 2;
    } else {
      copy_len = 14 + *len_stream++; // why is the value not 16 here, the above case copies up to 16 bytes.
      if ((uintptr_t)copy_len >(uintptr_t)(dst_end - dst))
        return false; // copy length out of bounds
      COPY_64(dst, copyfrom);
      COPY_64(dst + 8, copyfrom + 8);
      COPY_64(dst + 16, copyfrom + 16);
      do {
        COPY_64(dst + 24, copyfrom + 24);
        copy_len -= 8;
        dst += 8;
        copyfrom += 8;
      } while (copy_len > 24);
      dst += copy_len;
    }
  }

  // check for incorrect input
  if (distance_stream != distance_stream_end ||
      len_stream != litlen_stream_end)
    return false;

  final_len = dst_end - dst;
  if (final_len != lit_stream_end - lit_stream)
    return false;

  if (final_len >= 8) {
    do {
      COPY_64_ADD(dst, lit_stream, &dst[last_lit_dist]);
      dst += 8;
      lit_stream += 8;
      final_len -= 8;
    } while (final_len >= 8);
  }
  if (final_len > 0) {
    do {
      *dst = *lit_stream++ + dst[last_lit_dist];
      dst++;
    } while (--final_len);
  }
  return true;
}


// Note: may access memory out of bounds on invalid input.
bool Kraken_ProcessLzRuns_Type1(KrakenLzTable *lzt, byte *dst, byte *dst_end, byte *dst_start) {
  const byte *flag_stream = lzt->flag_stream, 
             *flag_stream_end = flag_stream + lzt->flag_stream_size;
  const int *len_stream = lzt->len_stream;
  const int *litlen_stream_end = lzt->len_stream + lzt->len_stream_size;
  const byte *lit_stream = lzt->lit_stream;
  const byte *lit_stream_end = lzt->lit_stream + lzt->lit_stream_size;
  const int *distance_stream = lzt->dist_stream;
  const int *distance_stream_end = lzt->dist_stream + lzt->dist_stream_size;
  const byte *copyfrom;
  uint32 final_len;
  int32 distance;
  int32 recent_distance[7];

  recent_distance[3] = -8;
  recent_distance[4] = -8;
  recent_distance[5] = -8;

  while (flag_stream < flag_stream_end) {
    uint32 f = *flag_stream++;
    uint32 litlen = f & 3;
    uint32 disttype = f >> 6;
    uint32 copy_len = (f >> 2) & 0xF;
  
    // use cmov
    uint32 next_long_length = *len_stream;
    const int *next_len_stream = len_stream + 1;

    len_stream = (litlen == 3) ? next_len_stream : len_stream; 
    litlen = (litlen == 3) ? next_long_length : litlen;
    recent_distance[6] = *distance_stream;

    COPY_64(dst, lit_stream);
    if (litlen > 8) {
      COPY_64(dst + 8, lit_stream + 8);
      if (litlen > 16) {
        COPY_64(dst + 16, lit_stream + 16);
        if (litlen > 24) {
          do {
            COPY_64(dst + 24, lit_stream + 24);
            litlen -= 8;
            dst += 8;
            lit_stream += 8;
          } while (litlen > 24);
        }
      }
    }
    dst += litlen;
    lit_stream += litlen;

    distance = recent_distance[disttype + 3];
    recent_distance[disttype + 3] = recent_distance[disttype + 2];
    recent_distance[disttype + 2] = recent_distance[disttype + 1];
    recent_distance[disttype + 1] = recent_distance[disttype + 0];
    recent_distance[3] = distance;
    
    distance_stream = (int*)((intptr_t)distance_stream + ((disttype + 1) & 4));

    if ((uintptr_t)distance < (uintptr_t)(dst_start - dst))
      return false; // distance out of bounds

    copyfrom = dst + distance;
    if (copy_len != 15) {
      COPY_64(dst, copyfrom);
      COPY_64(dst + 8, copyfrom + 8);
      dst += copy_len + 2;
    } else {
      copy_len = 14 + *len_stream++; // why is the value not 16 here, the above case copies up to 16 bytes.
      if ((uintptr_t)copy_len > (uintptr_t)(dst_end - dst))
        return false; // copy length out of bounds
      COPY_64(dst, copyfrom);
      COPY_64(dst + 8, copyfrom + 8);
      COPY_64(dst + 16, copyfrom + 16);
      do {
        COPY_64(dst + 24, copyfrom + 24);
        copy_len -= 8;
        dst += 8;
        copyfrom += 8;
      } while (copy_len > 24);
      dst += copy_len;
    }
  }

  // check for incorrect input
  if (distance_stream != distance_stream_end ||
      len_stream != litlen_stream_end)
    return false;

  final_len = dst_end - dst;
  if (final_len != lit_stream_end - lit_stream)
    return false;

  if (final_len >= 8) {
    do {
      COPY_64(dst, lit_stream);
      dst += 8;
      lit_stream += 8;
      final_len -= 8;
    } while (final_len >= 8);
  }
  if (final_len > 0) {
    do {
      *dst++ = *lit_stream++;
    } while (--final_len);
  }
  return true;
}

bool Kraken_ProcessLzRuns(int mode, byte *dst, int dst_size, int offset, KrakenLzTable *lztable) {
  byte *dst_end = dst + dst_size;

  if (mode == 1)
    return Kraken_ProcessLzRuns_Type1(lztable, dst + (offset == 0 ? 8 : 0), dst_end, dst - offset);

  if (mode == 0)
    return Kraken_ProcessLzRuns_Type0(lztable, dst + (offset == 0 ? 8 : 0), dst_end, dst - offset);


  return false;
}

// Decode one 256kb big quantum block. It's divided into two 128k blocks
// internally that are compressed separately but with a shared history.
int Kraken_DecodeQuantum(byte *dst, byte *dst_end, byte *dst_start,
                         const byte *src, const byte *src_end,
                         byte *temp, byte *temp_end) {
  const byte *src_in = src;
  int mode, chunkhdr, dst_count, src_used, written_bytes;

  while (dst_end - dst != 0) {
    dst_count = dst_end - dst;
    if (dst_count > 0x20000) dst_count = 0x20000;
    if (src_end - src < 4)
      return -1;
    chunkhdr = src[2] | src[1] << 8 | src[0] << 16;
    if (!(chunkhdr & 0x800000)) {
      // Stored as plain huffman without any match copying.
      byte *out = dst;
      src_used = Kraken_DecodeBytes(&out, src, src_end, &written_bytes, dst_count);
      if (src_used < 0 || written_bytes != dst_count)
        return -1;
    } else {
      src += 3;
      src_used = chunkhdr & 0x7FFFF;
      mode = (chunkhdr >> 19) & 0xF;
      if (src_end - src < src_used)
        return -1;
      if (src_used < dst_count) {
        int temp_usage = 2 * dst_count + 32;
        if (temp_usage > 0x40000) temp_usage = 0x40000;
        if (!Kraken_ReadLzTable(mode,
                               src, src + src_used,
                               dst, dst_count,
                               dst - dst_start,
                               temp + sizeof(KrakenLzTable), temp + temp_usage,
                               (KrakenLzTable*)temp))
          return -1;
        if (!Kraken_ProcessLzRuns(mode, dst, dst_count, dst - dst_start, (KrakenLzTable*)temp))
          return -1;
      } else if (src_used > dst_count || mode != 0) {
        return -1;
      } else {
        memmove(dst, src, dst_count);
      }
    }
    src += src_used;
    dst += dst_count;
  }
  return src - src_in;
}


int Mermaid_DecodeFarOffsets(const byte *src, const byte *src_end, uint32 *output, size_t output_size, int64 offset) {
  const byte *src_cur = src;
  size_t i;
  uint32 off;

  if (offset < (0xC00000 - 1)) {
    for (i = 0; i != output_size; i++) {
      if (src_end - src_cur < 3)
        return -1;
      off = src_cur[0] | src_cur[1] << 8 | src_cur[2] << 16;
      src_cur += 3;
      output[i] = off;
      if (off > offset)
        return -1;
    }
    return src_cur - src;
  }

  for (i = 0; i != output_size; i++) {
    if (src_end - src_cur < 3)
      return -1;
    off = src_cur[0] | src_cur[1] << 8 | src_cur[2] << 16;
    src_cur += 3;

    if (off >= 0xc00000) {
      if (src_cur == src_end)
        return -1;
      off += *src_cur++ << 22;
    }
    output[i] = off;
    if (off > offset)
      return -1;
  }
  return src_cur - src;
}


bool Mermaid_ReadLzTable(int mode,
                         const byte *src, const byte *src_end,
                         byte *dst, int dst_size, int64 offset,
                         byte *output, byte *output_end, MermaidLzTable *lz) {
  byte *out;
  int decode_count, n;
  uint32 tmp, off32_size_2, off32_size_1;

  if (mode > 1)
    return false;

  if (src_end - src < 10)
    return false;

  if (offset == 0) {
    COPY_64(dst, src);
    dst += 8;
    src += 8;
  }

  // Decode lit stream
  out = output;
  n = Kraken_DecodeBytes(&out, src, src_end, &decode_count, output_end - output);
  if (n < 0)
    return false;
  src += n;
  lz->lit_stream = out;
  lz->lit_stream_end = out + decode_count;
  output += decode_count;

  // Decode flag stream
  out = output;
  n = Kraken_DecodeBytes(&out, src, src_end, &decode_count, output_end - output);
  if (n < 0)
    return false;
  src += n;
  lz->flag_stream = out;
  lz->flag_stream_end = out + decode_count;
  output += decode_count;
  
  lz->flag_stream_2_offs_end = decode_count;
  if (dst_size <= 0x10000) {
    lz->flag_stream_2_offs = decode_count;
  } else {
    if (src_end - src < 2)
      return false;
    lz->flag_stream_2_offs = *(uint16*)src;
    src += 2;
  }

  lz->off16_stream = (uint16*)(src + 2);
  src += 2 + *(uint16*)src * 2;
  lz->off16_stream_end = (uint16*)src;

  if (src_end - src < 3)
    return false;

  tmp = src[0] | src[1] << 8 | src[2] << 16;
  src += 3;

  if (tmp != 0) {
    off32_size_1 = tmp >> 12;
    off32_size_2 = tmp & 0xFFF;
    if (off32_size_1 == 4095) {
      if (src_end - src < 2)
        return false;
      off32_size_1 = *(uint16*)src;
      src += 2;
    }
    if (off32_size_2 == 4095) {
      if (src_end - src < 2)
        return false;
      off32_size_2 = *(uint16*)src;
      src += 2;
    }
    lz->off32_size_1 = off32_size_1;
    lz->off32_size_2 = off32_size_2;

    if (output_end - output < 4 * (off32_size_2 + off32_size_1) + 64)
      return false;

    lz->off32_stream_1 = (uint32*)output;
    output += off32_size_1 * 4;
    // store dummy bytes after for prefetcher.
    ((uint64*)output)[0] = 0;
    ((uint64*)output)[1] = 0;
    ((uint64*)output)[2] = 0;
    ((uint64*)output)[3] = 0;
    output += 32;

    lz->off32_stream_2 = (uint32*)output;
    output += off32_size_2 * 4;
    // store dummy bytes after for prefetcher.
    ((uint64*)output)[0] = 0;
    ((uint64*)output)[1] = 0;
    ((uint64*)output)[2] = 0;
    ((uint64*)output)[3] = 0;
    output += 32;

    n = Mermaid_DecodeFarOffsets(src, src_end, lz->off32_stream_1, lz->off32_size_1, offset);
    if (n < 0)
      return false;
    src += n;

    n = Mermaid_DecodeFarOffsets(src, src_end, lz->off32_stream_2, lz->off32_size_2, offset + 0x10000);
    if (n < 0)
      return false;
    src += n;
  } else {
    if (output_end - output < 32)
      return false;
    lz->off32_size_1 = 0;
    lz->off32_size_2 = 0;
    lz->off32_stream_1 = (uint32*)output;
    lz->off32_stream_2 = (uint32*)output;
    // store dummy bytes after for prefetcher.
    ((uint64*)output)[0] = 0;
    ((uint64*)output)[1] = 0;
    ((uint64*)output)[2] = 0;
    ((uint64*)output)[3] = 0;
  }
  lz->length_stream = src;
  return true;
}

const byte *Mermaid_Mode0(byte *dst, size_t dst_size, byte *dst_ptr_end, byte *dst_start,
                          const byte *src_end, MermaidLzTable *lz, int32 *saved_dist, size_t startoff) {
  const byte *dst_end = dst + dst_size;
  const byte *flag_stream = lz->flag_stream;
  const byte *flag_stream_end = lz->flag_stream_end;
  const byte *length_stream = lz->length_stream;
  const byte *lit_stream = lz->lit_stream;
  const byte *lit_stream_end = lz->lit_stream_end;
  const uint16 *off16_stream = lz->off16_stream;
  const uint16 *off16_stream_end = lz->off16_stream_end;
  const uint32 *off32_stream = lz->off32_stream;
  const uint32 *off32_stream_end = lz->off32_stream_end;
  intptr_t recent_distance = *saved_dist;
  const byte *match;
  intptr_t length;
  const byte *dst_begin = dst;

  dst += startoff;

  while (flag_stream < flag_stream_end) {
    uintptr_t flag = *flag_stream++;
    if (flag >= 24) {
      intptr_t new_dist = *off16_stream;
      uintptr_t use_distance = (uintptr_t)(flag >> 7) - 1;
      uintptr_t litlen = (flag & 7);
      COPY_64_ADD(dst, lit_stream, &dst[recent_distance]);
      dst += litlen;
      lit_stream += litlen;
      recent_distance ^= use_distance & (recent_distance ^ -new_dist);
      off16_stream = (uint16*)((uintptr_t)off16_stream + (use_distance & 2));
      match = dst + recent_distance;
      COPY_64(dst, match);
      COPY_64(dst + 8, match + 8);
      dst += (flag >> 3) & 0xF;
    } else if (flag > 2) {
      length = flag + 5;

      if (off32_stream == off32_stream_end)
        return NULL;
      match = dst_begin - *off32_stream++;
      recent_distance = (match - dst);

      if (dst_end - dst < length)
        return NULL;
      COPY_64(dst, match);
      COPY_64(dst + 8, match + 8);
      COPY_64(dst + 16, match + 16);
      COPY_64(dst + 24, match + 24);
      dst += length;
      _mm_prefetch((char*)dst_begin - off32_stream[3], _MM_HINT_T0);
    } else if (flag == 0) {
      if (src_end - length_stream == 0)
        return NULL;
      length = *length_stream;
      if (length > 251) {
        if (src_end - length_stream < 3)
          return NULL;
        length += (size_t)*(uint16*)(length_stream + 1) * 4;
        length_stream += 2;
      }
      length_stream += 1;

      length += 64;
      if (dst_end - dst < length ||
          lit_stream_end - lit_stream < length)
        return NULL;

      do {
        COPY_64_ADD(dst, lit_stream, &dst[recent_distance]);
        COPY_64_ADD(dst + 8, lit_stream + 8, &dst[recent_distance + 8]);
        dst += 16;
        lit_stream += 16;
        length -= 16;
      } while (length > 0);
      dst += length;
      lit_stream += length;
    } else if (flag == 1) {
      if (src_end - length_stream == 0)
        return NULL;
      length = *length_stream;
      if (length > 251) {
        if (src_end - length_stream < 3)
          return NULL;
        length += (size_t)*(uint16*)(length_stream + 1) * 4;
        length_stream += 2;
      }
      length_stream += 1;
      length += 91;

      if (off16_stream == off16_stream_end)
        return NULL;
      match = dst - *off16_stream++;
      recent_distance = (match - dst);
      do {
        COPY_64(dst, match);
        COPY_64(dst + 8, match + 8);
        dst += 16;
        match += 16;
        length -= 16;
      } while (length > 0);
      dst += length;
    } else /* flag == 2 */ {
      if (src_end - length_stream == 0)
        return NULL;
      length = *length_stream;
      if (length > 251) {
        if (src_end - length_stream < 3)
          return NULL;
        length += (size_t)*(uint16*)(length_stream + 1) * 4;
        length_stream += 2;
      }
      length_stream += 1;
      length += 29;
      if (off32_stream == off32_stream_end)
        return NULL;
      match = dst_begin - *off32_stream++;
      recent_distance = (match - dst);
      do {
        COPY_64(dst, match);
        COPY_64(dst + 8, match + 8);
        dst += 16;
        match += 16;
        length -= 16;
      } while (length > 0);
      dst += length;
      _mm_prefetch((char*)dst_begin - off32_stream[3], _MM_HINT_T0);
    }
  }

  length = dst_end - dst;
  if (length >= 8) {
    do {
      COPY_64_ADD(dst, lit_stream, &dst[recent_distance]);
      dst += 8;
      lit_stream += 8;
      length -= 8;
    } while (length >= 8);
  }
  if (length > 0) {
    do {
      *dst = *lit_stream++ + dst[recent_distance];
      dst++;
    } while (--length);
  }

  *saved_dist = (int32)recent_distance;
  lz->length_stream = length_stream;
  lz->off16_stream = off16_stream;
  lz->lit_stream = lit_stream;
  return length_stream;
}

const byte *Mermaid_Mode1(byte *dst, size_t dst_size, byte *dst_ptr_end, byte *dst_start,
                         const byte *src_end, MermaidLzTable *lz, int32 *saved_dist, size_t startoff) {
  const byte *dst_end = dst + dst_size;
  const byte *flag_stream = lz->flag_stream;
  const byte *flag_stream_end = lz->flag_stream_end;
  const byte *length_stream = lz->length_stream;
  const byte *lit_stream = lz->lit_stream;
  const byte *lit_stream_end = lz->lit_stream_end;
  const uint16 *off16_stream = lz->off16_stream;
  const uint16 *off16_stream_end = lz->off16_stream_end;
  const uint32 *off32_stream = lz->off32_stream;
  const uint32 *off32_stream_end = lz->off32_stream_end;
  intptr_t recent_distance = *saved_dist;
  const byte *match;
  intptr_t length;
  const byte *dst_begin = dst;

  dst += startoff;

  while (flag_stream < flag_stream_end) {
    uintptr_t flag = *flag_stream++;
    if (flag >= 24) {
      intptr_t new_dist = *off16_stream;
      uintptr_t use_distance = (uintptr_t)(flag >> 7) - 1;
      uintptr_t litlen = (flag & 7);
      COPY_64(dst, lit_stream);
      dst += litlen;
      lit_stream += litlen;
      recent_distance ^= use_distance & (recent_distance ^ -new_dist);
      off16_stream = (uint16*)((uintptr_t)off16_stream + (use_distance & 2));
      match = dst + recent_distance;
      COPY_64(dst, match);
      COPY_64(dst + 8, match + 8);
      dst += (flag >> 3) & 0xF;
    } else if (flag > 2) {
      length = flag + 5;

      if (off32_stream == off32_stream_end)
        return NULL;
      match = dst_begin - *off32_stream++;
      recent_distance = (match - dst);
      
      if (dst_end - dst < length)
        return NULL;
      COPY_64(dst, match);
      COPY_64(dst + 8, match + 8);
      COPY_64(dst + 16, match + 16);
      COPY_64(dst + 24, match + 24);
      dst += length;
      _mm_prefetch((char*)dst_begin - off32_stream[3], _MM_HINT_T0);
    } else if (flag == 0) {
      if (src_end - length_stream == 0)
        return NULL;
      length = *length_stream;
      if (length > 251) {
        if (src_end - length_stream < 3)
          return NULL;
        length += (size_t)*(uint16*)(length_stream + 1) * 4;
        length_stream += 2;
      }
      length_stream += 1;

      length += 64;
      if (dst_end - dst < length ||
          lit_stream_end - lit_stream < length)
        return NULL;

      do {
        COPY_64(dst, lit_stream);
        COPY_64(dst + 8, lit_stream + 8);
        dst += 16;
        lit_stream += 16;
        length -= 16;
      } while (length > 0);
      dst += length;
      lit_stream += length;
    } else if (flag == 1) {
      if (src_end - length_stream == 0)
        return NULL;
      length = *length_stream;
      if (length > 251) {
        if (src_end - length_stream < 3)
          return NULL;
        length += (size_t)*(uint16*)(length_stream + 1) * 4;
        length_stream += 2;
      }
      length_stream += 1;
      length += 91;
      
      if (off16_stream == off16_stream_end)
        return NULL;
      match = dst - *off16_stream++;
      recent_distance = (match - dst);
      do {
        COPY_64(dst, match);
        COPY_64(dst + 8, match + 8);
        dst += 16;
        match += 16;
        length -= 16;
      } while (length > 0);
      dst += length;
    } else /* flag == 2 */ {
      if (src_end - length_stream == 0)
        return NULL;
      length = *length_stream;
      if (length > 251) {
        if (src_end - length_stream < 3)
          return NULL;
        length += (size_t)*(uint16*)(length_stream + 1) * 4;
        length_stream += 2;
      }
      length_stream += 1;
      length += 29;

      if (off32_stream == off32_stream_end)
        return NULL;
      match = dst_begin - *off32_stream++;
      recent_distance = (match - dst);
      
      do {
        COPY_64(dst, match);
        COPY_64(dst + 8, match + 8);
        dst += 16;
        match += 16;
        length -= 16;
      } while (length > 0);
      dst += length;

      _mm_prefetch((char*)dst_begin - off32_stream[3], _MM_HINT_T0);
    }
  }

  length = dst_end - dst;
  if (length >= 8) {
    do {
      COPY_64(dst, lit_stream);
      dst += 8;
      lit_stream += 8;
      length -= 8;
    } while (length >= 8);
  }
  if (length > 0) {
    do {
      *dst++ = *lit_stream++;
    } while (--length);
  }

  *saved_dist = (int32)recent_distance;
  lz->length_stream = length_stream;
  lz->off16_stream = off16_stream;
  lz->lit_stream = lit_stream;
  return length_stream;
}

bool Mermaid_ProcessLzRuns(int mode,
                           const byte *src, const byte *src_end,
                           byte *dst, size_t dst_size, uint64 offset, byte *dst_end,
                           MermaidLzTable *lz) {
  
  int iteration = 0;
  byte *dst_start = dst - offset;
  int32 saved_dist = -8;
  const byte *src_cur;

  for (iteration = 0; iteration != 2; iteration++) {
    size_t dst_size_cur = dst_size;
    if (dst_size_cur > 0x10000) dst_size_cur = 0x10000;
    
    if (iteration == 0) {
      lz->off32_stream = lz->off32_stream_1;
      lz->off32_stream_end = lz->off32_stream_1 + lz->off32_size_1 * 4;
      lz->flag_stream_end = lz->flag_stream + lz->flag_stream_2_offs;
    } else {
      lz->off32_stream = lz->off32_stream_2;
      lz->off32_stream_end = lz->off32_stream_2 + lz->off32_size_2 * 4;
      lz->flag_stream_end = lz->flag_stream + lz->flag_stream_2_offs_end;
      lz->flag_stream += lz->flag_stream_2_offs;
    }

    if (mode == 0) {
      src_cur = Mermaid_Mode0(dst, dst_size_cur, dst_end, dst_start, src_end, lz, &saved_dist, 
        (offset == 0) && (iteration == 0) ? 8 : 0);
    } else {
      src_cur = Mermaid_Mode1(dst, dst_size_cur, dst_end, dst_start, src_end, lz, &saved_dist,
        (offset == 0) && (iteration == 0) ? 8 : 0);
    }
    if (src_cur == NULL)
      return false;

    dst += dst_size_cur;
    dst_size -= dst_size_cur;
    if (dst_size == 0)
      break;
  }

  if (src_cur != src_end)
    return false;


  return true;
}


int Mermaid_DecodeQuantum(byte *dst, byte *dst_end, byte *dst_start,
                          const byte *src, const byte *src_end,
                          byte *temp, byte *temp_end) {

  const byte *src_in = src;
  int mode, chunkhdr, dst_count, src_used, written_bytes;

  while (dst_end - dst != 0) {
    dst_count = dst_end - dst;
    if (dst_count > 0x20000) dst_count = 0x20000;
    if (src_end - src < 4)
      return -1;
    chunkhdr = src[2] | src[1] << 8 | src[0] << 16;
    if (!(chunkhdr & 0x800000)) {
      // Stored without any match copying.
      byte *out = dst;
      src_used = Kraken_DecodeBytes(&out, src, src_end, &written_bytes, dst_count);
      if (src_used < 0 || written_bytes != dst_count)
        return -1;
    } else {
      src += 3;
      src_used = chunkhdr & 0x7FFFF;
      mode = (chunkhdr >> 19) & 0xF;
      if (src_end - src < src_used)
        return -1;
      if (src_used < dst_count) {
        int temp_usage = 2 * dst_count + 32;
        if (temp_usage > 0x40000) temp_usage = 0x40000;
        if (!Mermaid_ReadLzTable(mode,
                                src, src + src_used,
                                dst, dst_count,
                                dst - dst_start,
                                temp + sizeof(MermaidLzTable), temp + temp_usage,
                                (MermaidLzTable*)temp))
          return -1;
        if (!Mermaid_ProcessLzRuns(mode,
                                   src, src + src_used,
                                   dst, dst_count,
                                   dst - dst_start, dst_end,
                                   (MermaidLzTable*)temp))
          return -1;
      } else if (src_used > dst_count || mode != 0) {
        return -1;
      } else {
        memmove(dst, src, dst_count);
      }
    }
    src += src_used;
    dst += dst_count;
  }
  return src - src_in;
}

bool Kraken_DecodeStep(struct KrakenDecoder *dec,
                       byte *dst_start, int offset, int dst_bytes_left,
                       const byte *src, int src_bytes_left) {
  const byte *src_in = src;
  const byte *src_end = src + src_bytes_left;
  KrakenQuantumHeader qhdr;
  KrakenHeader hdr;
  int n;

  src = Kraken_ParseHeader(&hdr, src);
  if (!src)
    return false;

  if (dst_bytes_left > 0x40000)
    dst_bytes_left = 0x40000;

  if (hdr.uncompressed) {
    if (src_end - src < dst_bytes_left) {
      dec->src_used = dec->dst_used = 0;
      return true;
    }
    memmove(dst_start + offset, src, dst_bytes_left);
    dec->src_used = (src - src_in) + dst_bytes_left;
    dec->dst_used = dst_bytes_left;
    return true;
  }

  src = Kraken_ParseQuantumHeader(&qhdr, src, hdr.use_checksums);
  if (!src || src > src_end)
    return false;

  if ((uintptr_t)(src_end - src) < qhdr.compressed_size) {
    dec->src_used = dec->dst_used = 0;
    return true;
  }
  
  if (qhdr.compressed_size > (uint32)dst_bytes_left)
    return false;

  if (qhdr.compressed_size == 0) {
    memset(dst_start + offset, qhdr.checksum, dst_bytes_left);
    dec->src_used = (src - src_in);
    dec->dst_used = dst_bytes_left;
    return true;
  }

  if (hdr.use_checksums && 
     (Kraken_GetCrc(src, qhdr.compressed_size) & 0xFFFFFF) != qhdr.checksum)
    return false;

  if (qhdr.compressed_size == dst_bytes_left) {
    memmove(dst_start + offset, src, dst_bytes_left);
    dec->src_used = (src - src_in) + dst_bytes_left;
    dec->dst_used = dst_bytes_left;
    return true;
  }

  if (hdr.decoder_type == 6) {
    n = Kraken_DecodeQuantum(dst_start + offset, dst_start + offset + dst_bytes_left, dst_start,
                         src, src + qhdr.compressed_size,
                         dec->phase_buf, dec->phase_buf + dec->phase_buf_size);
  } else {
    n = Mermaid_DecodeQuantum(dst_start + offset, dst_start + offset + dst_bytes_left, dst_start,
                              src, src + qhdr.compressed_size,
                              dec->phase_buf, dec->phase_buf + dec->phase_buf_size);
  }

  if (n != qhdr.compressed_size)
    return false;

  dec->src_used = (src - src_in) + n;
  dec->dst_used = dst_bytes_left;
  return true;
}
  
int Kraken_Decompress(const byte *src, size_t src_len, byte *dst, size_t dst_len) {
  KrakenDecoder *dec = Kraken_Create();
  int offset = 0;
  while (dst_len != 0) {
    if (!Kraken_DecodeStep(dec, dst, offset, dst_len, src, src_len))
      goto FAIL;
    if (dec->src_used == 0)
      goto FAIL;
    src += dec->src_used;
    src_len -= dec->src_used;
    dst_len -= dec->dst_used;
    offset += dec->dst_used;
  }
  if (src_len != 0)
    goto FAIL;
  Kraken_Destroy(dec);
  return offset;
FAIL:
  Kraken_Destroy(dec);
  return -1;
}

// The decompressor will write outside of the target buffer.
#define SAFE_SPACE 64

void error(const char *s) {
  fprintf(stderr, "%s\n", s);
  exit(1);
}


byte *load_file(const char *filename, int *size) {
  FILE *f = fopen(filename, "rb");
  if (!f) error("file open error");
  fseek(f, 0, SEEK_END);
  int packed_size = ftell(f);
  fseek(f, 0, SEEK_SET);
  byte *input = new byte[packed_size];
  if (!input) error("memory error");
  if (fread(input, 1, packed_size, f) != packed_size) error("error reading");
  fclose(f);
  *size = packed_size;
  return input;
}


int main(int argc, char *argv[]) {
  __int64 start, end, freq;

  if (argc != 3) {
    fprintf(stderr, "unkraken v0.01\n");
    error("unkraken input output");
  }

  int packed_size;
  byte *input = load_file(argv[1], &packed_size);
  
  uint32 unpacked_size = *(uint32*)input;
  if (unpacked_size > 1024 * 1024 * 1024) error("file too large");

  byte *output = new byte[unpacked_size + SAFE_SPACE];
  if (!output) error("memory error");

  int outbytes;

  for (int i = 0; i < 1; i++) {
    QueryPerformanceCounter((LARGE_INTEGER*)&start);

    outbytes = Kraken_Decompress(input + 4, packed_size - 4, output, unpacked_size);
    if (outbytes != unpacked_size) error("size error");

    QueryPerformanceCounter((LARGE_INTEGER*)&end);
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

    double seconds = (double)(end - start) / freq;
    fprintf(stderr, "%-20s: %8d => %8d (%.2f seconds, %.2f MB/s)\n", argv[1], packed_size, unpacked_size, seconds, unpacked_size * 1e-6 / seconds);

  }

#if 0
  int ref_size;
  byte *ref = load_file(argv[2], &ref_size);

  if (ref_size != unpacked_size) error("ref size error");

  for (int i = 0; i != ref_size; i++) {
    if (ref[i] != output[i]) {
      fprintf(stderr, "Decode error at 0x%x. Was %d wanted %d\n", i, output[i], ref[i]);
      break;
    }
  }
#else
  FILE *f = fopen(argv[2], "wb");
  if (!f) error("file open for write error");
  fwrite(output, 1, outbytes, f);
  fclose(f);
#endif

  return 0;
}

