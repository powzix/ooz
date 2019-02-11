# ooz
Open source Kraken / Mermaid / Selkie / Leviathan / LZNA / Bitknit decompressor

Also supports using oo2core_7_win32.dll / oo2core_7_win64.dll which can be acquired from the free game Warframe on steam.

```
ooz v7.0

Usage: ooz [options] input [output]
 -c --stdout              write to stdout
 -d --decompress          decompress (default)
 -z --compress            compress (requires oo2core_7_win64.dll)
 -b                       just benchmark, don't overwrite anything
 -f                       force overwrite existing file
 --dll                    decompress with the dll
 --verify                 decompress and verify that it matches output
 --verify=<folder>        verify with files in this folder
 -<1-9> --level=<-4..10>  compression level
 -m<k>                    [k|m|s|l|h] compressor selection
 --kraken --mermaid --selkie --leviathan --hydra    compressor selection

(Warning! not fuzz safe, so please trust the input)
```
