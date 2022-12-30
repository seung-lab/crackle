Note that all numbers provided are for the "flat" label layout. The pin layout is still under investigation (it is much slower and may not always provide a benefit).

## Size Benchmark
shape: (256, 256, 64)
PINKY40 CUTOUTS (connectomics.npy)

      ckl:         40505    
      cpso:       104686  (38.69%) 
      raw:      16777216  (0.24%)
      ckl.gz:      29429   
      cpso.gz      42991  (68.45%)
      raw:.gz     103401  (28.46%)
    

      ckl:         11803    
      cpso:        23298  (50.66%) 
      raw:      16777216  (0.07%)
      ckl.gz:       5758   
      cpso.gz       6772  (85.03%)
      raw:.gz      23622  (24.38%)
    

      ckl:         21609    
      cpso:        49196  (43.92%) 
      raw:      16777216  (0.13%)
      ckl.gz:      14039   
      cpso.gz      19212  (73.07%)
      raw:.gz      56435  (24.88%)
    
WATERSHED CUTOUTS (ws.npy)

      ckl:        574187    
      cpso:      2758348  (20.82%) 
      raw:      33554432  (1.71%)
      ckl.gz:     431338   
      cpso.gz     703903  (61.28%)
      raw:.gz    1038695  (41.53%)
    

      ckl:        542143    
      cpso:      2534896  (21.39%) 
      raw:      33554432  (1.62%)
      ckl.gz:     404924   
      cpso.gz     656077  (61.72%)
      raw:.gz     984766  (41.12%)
    

      ckl:        564529    
      cpso:      2685326  (21.02%) 
      raw:      33554432  (1.68%)
      ckl.gz:     423111   
      cpso.gz     694604  (60.91%)
      raw:.gz    1010223  (41.88%)
    
RANDOM NOISE [0,2000) uint32

      ckl:       8405607    
      cpso:     17302002  (48.58%) 
      raw:      16777216  (50.10%)
      ckl.gz:    6792130   
      cpso.gz    8024626  (84.64%)
      raw:.gz    8023814  (84.65%)
    

      ckl:       8405644    
      cpso:     17302000  (48.58%) 
      raw:      16777216  (50.10%)
      ckl.gz:    6792245   
      cpso.gz    8024848  (84.64%)
      raw:.gz    8023877  (84.65%)
    

      ckl:       8406168    
      cpso:     17301992  (48.58%) 
      raw:      16777216  (50.10%)
      ckl.gz:    6791867   
      cpso.gz    8024867  (84.64%)
      raw:.gz    8023793  (84.65%)
    
BINARY NOISE [0,1] uint8 (pathological case)

      ckl:       1843818    
      cpso:      2900459  (63.57%) 
      raw:       4194304  (43.96%)
      ckl.gz:    1194228   
      cpso.gz    1098294  (108.73%)
      raw:.gz     667171  (179.00%)
    

      ckl:       1842842    
      cpso:      2898730  (63.57%) 
      raw:       4194304  (43.94%)
      ckl.gz:    1194008   
      cpso.gz    1099130  (108.63%)
      raw:.gz     667111  (178.98%)
    

      ckl:       1844494    
      cpso:      2899730  (63.61%) 
      raw:       4194304  (43.98%)
      ckl.gz:    1194865   
      cpso.gz    1100276  (108.60%)
      raw:.gz     666783  (179.20%)

EMPTY

      ckl:           607    
      cpso:          824  (73.67%) 
      raw:      16777216  (0.00%)
      ckl.gz:         34   
      cpso.gz         43  (79.07%)
      raw:.gz      16316  (0.21%)
    
SOLID ONES

      ckl:           607    
      cpso:          824  (73.67%) 
      raw:      16777216  (0.00%)
      ckl.gz:         36   
      cpso.gz         47  (76.60%)
      raw:.gz      16322  (0.22%)

## Performance Benchmark

shape: (256, 256, 64)  
PINKY40 CUTOUTS (connectomics.npy)

      compress     :  51.89 MVx/sec (206730 bytes, 1.2%)
      compress+gz  :  49.86 MVx/sec (139201 bytes, 67.3%)
      decompress   :  61.87 MVx/sec
      decompress+gz:  61.19 MVx/sec
    

      compress     :  90.83 MVx/sec (90864 bytes, 0.5%)
      compress+gz  :  88.39 MVx/sec (67940 bytes, 74.8%)
      decompress   :  72.77 MVx/sec
      decompress+gz:  72.38 MVx/sec
    

      compress     :  64.93 MVx/sec (143477 bytes, 0.9%)
      compress+gz  :  63.00 MVx/sec (109468 bytes, 76.3%)
      decompress   :  67.88 MVx/sec
      decompress+gz:  67.34 MVx/sec
    
WATERSHED CUTOUTS (ws.npy)

      compress     :  24.21 MVx/sec (573341 bytes, 1.7%)
      compress+gz  :  23.25 MVx/sec (428509 bytes, 74.7%)
      decompress   :  45.24 MVx/sec
      decompress+gz:  44.34 MVx/sec
    

      compress     :  25.16 MVx/sec (552975 bytes, 1.6%)
      compress+gz  :  24.17 MVx/sec (410592 bytes, 74.3%)
      decompress   :  49.42 MVx/sec
      decompress+gz:  48.43 MVx/sec
    

      compress     :  25.08 MVx/sec (556677 bytes, 1.7%)
      compress+gz  :  24.10 MVx/sec (414916 bytes, 74.5%)
      decompress   :  48.29 MVx/sec
      decompress+gz:  47.33 MVx/sec
    
RANDOM NOISE [0,2000) uint32

      compress     :  11.37 MVx/sec (8405666 bytes, 50.1%)
      compress+gz  :  9.03 MVx/sec (6882837 bytes, 81.9%)
      decompress   :  59.66 MVx/sec
      decompress+gz:  41.83 MVx/sec
    

      compress     :  11.31 MVx/sec (8405540 bytes, 50.1%)
      compress+gz  :  9.00 MVx/sec (6882446 bytes, 81.9%)
      decompress   :  70.85 MVx/sec
      decompress+gz:  47.18 MVx/sec
    

      compress     :  11.43 MVx/sec (8405455 bytes, 50.1%)
      compress+gz  :  9.07 MVx/sec (6882718 bytes, 81.9%)
      decompress   :  69.72 MVx/sec
      decompress+gz:  46.73 MVx/sec
    
BINARY NOISE [0,1] uint8 (pathological case)

      compress     :  6.26 MVx/sec (1841889 bytes, 43.9%)
      compress+gz  :  6.01 MVx/sec (1193359 bytes, 64.8%)
      decompress   :  21.98 MVx/sec
      decompress+gz:  21.43 MVx/sec
    

      compress     :  6.28 MVx/sec (1844669 bytes, 44.0%)
      compress+gz  :  6.03 MVx/sec (1194496 bytes, 64.8%)
      decompress   :  22.61 MVx/sec
      decompress+gz:  22.02 MVx/sec
    

      compress     :  6.27 MVx/sec (1843716 bytes, 44.0%)
      compress+gz  :  6.02 MVx/sec (1194203 bytes, 64.8%)
      decompress   :  22.07 MVx/sec
      decompress+gz:  21.51 MVx/sec
    
EMPTY

      compress     :  31.13 MVx/sec (607 bytes, 0.0%)
      compress+gz  :  31.12 MVx/sec (46 bytes, 7.6%)
      decompress   :  78.12 MVx/sec
      decompress+gz:  78.11 MVx/sec
    
SOLID ONES

      compress     :  29.83 MVx/sec (607 bytes, 0.0%)
      compress+gz  :  29.83 MVx/sec (48 bytes, 7.9%)
      decompress   :  82.59 MVx/sec
      decompress+gz:  82.58 MVx/sec
