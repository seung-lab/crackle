Note that all numbers provided are for the "flat" label layout. The pin layout is still under investigation (it is much slower and may not always provide a benefit). All data were first converted to Fortran order.

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

      compress     :  105.91 MVx/sec (58767 bytes, 0.4%)
      compress+gz  :  103.73 MVx/sec (43407 bytes, 73.9%)
      decompress   :  72.90 MVx/sec
      decompress+gz:  72.47 MVx/sec
    

      compress     :  57.63 MVx/sec (173198 bytes, 1.0%)
      compress+gz  :  55.79 MVx/sec (131285 bytes, 75.8%)
      decompress   :  66.40 MVx/sec
      decompress+gz:  65.79 MVx/sec
    

      compress     :  54.05 MVx/sec (219674 bytes, 1.3%)
      compress+gz  :  51.94 MVx/sec (148645 bytes, 67.7%)
      decompress   :  64.94 MVx/sec
      decompress+gz:  64.29 MVx/sec
    
WATERSHED CUTOUTS (ws.npy)

      compress     :  24.13 MVx/sec (585525 bytes, 1.7%)
      compress+gz  :  23.18 MVx/sec (436503 bytes, 74.5%)
      decompress   :  47.94 MVx/sec
      decompress+gz:  46.96 MVx/sec
    

      compress     :  24.89 MVx/sec (558694 bytes, 1.7%)
      compress+gz  :  23.91 MVx/sec (418787 bytes, 75.0%)
      decompress   :  48.52 MVx/sec
      decompress+gz:  47.55 MVx/sec
    

      compress     :  24.88 MVx/sec (562736 bytes, 1.7%)
      compress+gz  :  23.90 MVx/sec (419293 bytes, 74.5%)
      decompress   :  47.31 MVx/sec
      decompress+gz:  46.35 MVx/sec
    
RANDOM NOISE [0,2000) uint32

      compress     :  14.65 MVx/sec (8405711 bytes, 50.1%)
      compress+gz  :  10.98 MVx/sec (6882592 bytes, 81.9%)
      decompress   :  60.78 MVx/sec
      decompress+gz:  42.54 MVx/sec
    

      compress     :  15.02 MVx/sec (8405429 bytes, 50.1%)
      compress+gz  :  11.18 MVx/sec (6882826 bytes, 81.9%)
      decompress   :  68.45 MVx/sec
      decompress+gz:  46.16 MVx/sec
    

      compress     :  14.78 MVx/sec (8405280 bytes, 50.1%)
      compress+gz  :  11.05 MVx/sec (6882707 bytes, 81.9%)
      decompress   :  66.84 MVx/sec
      decompress+gz:  45.45 MVx/sec
    
BINARY NOISE [0,1] uint8 (pathological case)

      compress     :  6.57 MVx/sec (1840310 bytes, 43.9%)
      compress+gz  :  6.29 MVx/sec (1191938 bytes, 64.8%)
      decompress   :  20.60 MVx/sec
      decompress+gz:  20.11 MVx/sec
    

      compress     :  6.61 MVx/sec (1840982 bytes, 43.9%)
      compress+gz  :  6.34 MVx/sec (1193078 bytes, 64.8%)
      decompress   :  21.59 MVx/sec
      decompress+gz:  21.05 MVx/sec
    

      compress     :  6.60 MVx/sec (1841117 bytes, 43.9%)
      compress+gz  :  6.32 MVx/sec (1192435 bytes, 64.8%)
      decompress   :  21.35 MVx/sec
      decompress+gz:  20.82 MVx/sec
    
EMPTY

      compress     :  170.79 MVx/sec (607 bytes, 0.0%)
      compress+gz  :  170.62 MVx/sec (46 bytes, 7.6%)
      decompress   :  74.61 MVx/sec
      decompress+gz:  74.60 MVx/sec
    
SOLID ONES

      compress     :  210.25 MVx/sec (607 bytes, 0.0%)
      compress+gz  :  210.02 MVx/sec (48 bytes, 7.9%)
      decompress   :  81.08 MVx/sec
      decompress+gz:  81.07 MVx/sec
