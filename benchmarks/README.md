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

      compress     :  60.78 MVx/sec (91880 bytes, 0.5%)
      compress+gz  :  59.67 MVx/sec (70870 bytes, 77.1%)
      decompress   :  67.10 MVx/sec
      decompress+gz:  66.59 MVx/sec
    

      compress     :  48.26 MVx/sec (137213 bytes, 0.8%)
      compress+gz  :  47.23 MVx/sec (102603 bytes, 74.8%)
      decompress   :  67.13 MVx/sec
      decompress+gz:  66.62 MVx/sec
    

      compress     :  90.90 MVx/sec (54905 bytes, 0.3%)
      compress+gz  :  89.48 MVx/sec (40614 bytes, 74.0%)
      decompress   :  75.05 MVx/sec
      decompress+gz:  74.80 MVx/sec
    
WATERSHED CUTOUTS (ws.npy)

      compress     :  16.07 MVx/sec (562503 bytes, 1.7%)
      compress+gz  :  15.65 MVx/sec (420223 bytes, 74.7%)
      decompress   :  46.31 MVx/sec
      decompress+gz:  45.42 MVx/sec
    

      compress     :  16.28 MVx/sec (559439 bytes, 1.7%)
      compress+gz  :  15.86 MVx/sec (416760 bytes, 74.5%)
      decompress   :  47.73 MVx/sec
      decompress+gz:  46.79 MVx/sec
    

      compress     :  16.30 MVx/sec (555993 bytes, 1.7%)
      compress+gz  :  15.88 MVx/sec (416685 bytes, 74.9%)
      decompress   :  46.92 MVx/sec
      decompress+gz:  46.01 MVx/sec
    
RANDOM NOISE [0,2000) uint32

      compress     :  9.24 MVx/sec (8405434 bytes, 50.1%)
      compress+gz  :  7.64 MVx/sec (6882581 bytes, 81.9%)
      decompress   :  64.59 MVx/sec
      decompress+gz:  43.61 MVx/sec
    

      compress     :  9.48 MVx/sec (8405648 bytes, 50.1%)
      compress+gz  :  7.80 MVx/sec (6882623 bytes, 81.9%)
      decompress   :  70.60 MVx/sec
      decompress+gz:  47.14 MVx/sec
    

      compress     :  9.54 MVx/sec (8405767 bytes, 50.1%)
      compress+gz  :  7.84 MVx/sec (6882553 bytes, 81.9%)
      decompress   :  70.71 MVx/sec
      decompress+gz:  47.17 MVx/sec
    

BINARY NOISE [0,1] uint8 (pathological case)

      compress     :  4.31 MVx/sec (1843317 bytes, 43.9%)
      compress+gz  :  4.19 MVx/sec (1193619 bytes, 64.8%)
      decompress   :  21.72 MVx/sec
      decompress+gz:  21.17 MVx/sec
    

      compress     :  4.33 MVx/sec (1844347 bytes, 44.0%)
      compress+gz  :  4.21 MVx/sec (1194077 bytes, 64.7%)
      decompress   :  21.51 MVx/sec
      decompress+gz:  20.97 MVx/sec
    

      compress     :  4.33 MVx/sec (1842973 bytes, 43.9%)
      compress+gz  :  4.21 MVx/sec (1193080 bytes, 64.7%)
      decompress   :  21.50 MVx/sec
      decompress+gz:  20.96 MVx/sec

EMPTY

      compress     :  26.56 MVx/sec (607 bytes, 0.0%)
      compress+gz  :  26.55 MVx/sec (46 bytes, 7.6%)
      decompress   :  82.14 MVx/sec
      decompress+gz:  82.13 MVx/sec
    
SOLID ONES

      compress     :  35.59 MVx/sec (607 bytes, 0.0%)
      compress+gz  :  35.59 MVx/sec (48 bytes, 7.9%)
      decompress   :  82.36 MVx/sec
      decompress+gz:  82.35 MVx/sec

