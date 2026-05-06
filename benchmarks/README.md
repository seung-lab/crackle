Note that all numbers provided are for the "flat" label layout. The pin layout is still under investigation (it is much slower and may not always provide a benefit). All data were first converted to Fortran order.

Machine: Mabook Pro M3

## Size Benchmark
shape: (256, 256, 64)  

### PINKY40 CUTOUTS (connectomics.npy)

      ckl:          94049   (0.56%)
      ckl(pin):     85699   (0.51%) 
      ckl(mkv):     76198   (0.45%)
      cpso:        258208  (1.54%)
      raw:       16777216  (100.00%)

      ckl.gz:          71828   (27.47%)
      ckl.gz(pin):     68383   (26.15%) 
      ckl.gz(mkv):     66263   (25.34%)
      cpso.gz:        110145  (42.12%)
      raw.gz:         261496  (100.00%)

      ckl.br:          66497   (41.54%)
      ckl.br(pin):     63759   (39.83%) 
      ckl.br(mkv):     63886   (39.91%)
      cpso.br:         92544  (57.81%)
      raw.br:         160076  (100.00%)

### BINARY IMAGE (connectomics.npy, label 67699431)

      ckl:          43401   (0.03%)
      ckl(pin):     41996   (0.03%) 
      ckl(mkv):     34886   (0.03%)
      cpso:         87480  (0.07%)
      raw:       134217728  (100.00%)

      ckl.gz:          24070   (11.62%)
      ckl.gz(pin):     23862   (11.52%) 
      ckl.gz(mkv):     22577   (10.90%)
      cpso.gz:         45130  (21.78%)
      raw.gz:         207188  (100.00%)

      ckl.br:          21552   (29.25%)
      ckl.br(pin):     21434   (29.09%) 
      ckl.br(mkv):     21094   (28.63%)
      cpso.br:         39110  (53.09%)
      raw.br:          73672  (100.00%)

### WATERSHED CUTOUTS (ws.npy)

      ckl:         553454   (1.65%)
      ckl(pin):    543172   (1.62%) 
      ckl(mkv):    492454   (1.47%)
      cpso:       2546412  (7.59%)
      raw:       33554432  (100.00%)

      ckl.gz:         409447   (40.89%)
      ckl.gz(pin):    403366   (40.28%) 
      ckl.gz(mkv):    396406   (39.58%)
      cpso.gz:        657185  (65.62%)
      raw.gz:        1001437  (100.00%)

      ckl.br:         370373   (68.99%)
      ckl.br(pin):    370793   (69.07%) 
      ckl.br(mkv):    369504   (68.82%)
      cpso.br:        543861  (101.30%)
      raw.br:         536875  (100.00%)

### EMPTY

      ckl:            993   (0.01%)
      ckl(pin):       930   (0.01%) 
      ckl(mkv):       993   (0.01%)
      cpso:           824  (0.00%)
      raw:       16777216  (100.00%)

      ckl.gz:             44   (0.27%)
      ckl.gz(pin):        42   (0.26%) 
      ckl.gz(mkv):        44   (0.27%)
      cpso.gz:            43  (0.26%)
      raw.gz:          16316  (100.00%)

      ckl.br:             40   (148.15%)
      ckl.br(pin):        40   (148.15%) 
      ckl.br(mkv):        40   (148.15%)
      cpso.br:            46  (170.37%)
      raw.br:             27  (100.00%)

### SOLID ONES

      ckl:            993   (0.01%)
      ckl(pin):       930   (0.01%) 
      ckl(mkv):       993   (0.01%)
      cpso:           824  (0.00%)
      raw:       16777216  (100.00%)

      ckl.gz:             46   (0.28%)
      ckl.gz(pin):        44   (0.27%) 
      ckl.gz(mkv):        46   (0.28%)
      cpso.gz:            47  (0.29%)
      raw.gz:          16322  (100.00%)

      ckl.br:             41   (141.38%)
      ckl.br(pin):        39   (134.48%) 
      ckl.br(mkv):        41   (141.38%)
      cpso.br:            54  (186.21%)
      raw.br:             29  (100.00%)

### RANDOM NOISE [0,2000) uint32

      ckl:        8412099   (50.14%)
      ckl(pin):   8412099   (50.14%) 
      ckl(mkv):   8411377   (50.14%)
      cpso:      17302012  (103.13%)
      raw:       16777216  (100.00%)

      ckl.gz:        6791243   (84.64%)
      ckl.gz(pin):   6791243   (84.64%) 
      ckl.gz(mkv):   6791246   (84.64%)
      cpso.gz:       8024331  (100.01%)
      raw.gz:        8023874  (100.00%)

      ckl.br:        5830287   (90.63%)
      ckl.br(pin):   5830287   (90.63%) 
      ckl.br(mkv):   5829830   (90.62%)
      cpso.br:       6571752  (102.16%)
      raw.br:        6433114  (100.00%)

### BINARY NOISE [0,1] uint8 (pathological case)

      ckl:         284350   (43.39%)
      ckl(pin):    284350   (43.39%) 
      ckl(mkv):    270090   (41.21%)
      cpso:        468337  (71.46%)
      raw:         655360  (100.00%)

      ckl.gz:         181347   (173.94%)
      ckl.gz(pin):    181347   (173.94%) 
      ckl.gz(mkv):    181795   (174.37%)
      cpso.gz:        186901  (179.27%)
      raw.gz:         104259  (100.00%)

      ckl.br:         170516   (208.09%)
      ckl.br(pin):    170516   (208.09%) 
      ckl.br(mkv):    172139   (210.07%)
      cpso.br:        162174  (197.91%)
      raw.br:          81942  (100.00%)

## Performance Benchmark

shape: (256, 256, 64)
parallel: 1

PINKY40 CUTOUTS (connectomics.npy)

      compress     :  242.30 MVx/sec (118868 bytes, 0.7%)
      compress+gz  :  224.40 MVx/sec (91111 bytes, 76.6%)
      decompress   :  174.38 MVx/sec
      decompress+gz:  172.73 MVx/sec
    

      compress     :  402.26 MVx/sec (20854 bytes, 0.1%)
      compress+gz  :  394.13 MVx/sec (12304 bytes, 59.0%)
      decompress   :  198.22 MVx/sec
      decompress+gz:  197.89 MVx/sec
    

      compress     :  227.40 MVx/sec (137089 bytes, 0.8%)
      compress+gz  :  210.01 MVx/sec (104692 bytes, 76.4%)
      decompress   :  169.80 MVx/sec
      decompress+gz:  168.09 MVx/sec
    
WATERSHED CUTOUTS (ws.npy)

      compress     :  90.59 MVx/sec (562800 bytes, 1.7%)
      compress+gz  :  81.46 MVx/sec (415300 bytes, 73.8%)
      decompress   :  121.73 MVx/sec
      decompress+gz:  118.33 MVx/sec
    

      compress     :  85.36 MVx/sec (574709 bytes, 1.7%)
      compress+gz  :  76.97 MVx/sec (426732 bytes, 74.3%)
      decompress   :  116.64 MVx/sec
      decompress+gz:  113.49 MVx/sec
    

      compress     :  95.73 MVx/sec (553882 bytes, 1.7%)
      compress+gz  :  85.85 MVx/sec (409469 bytes, 73.9%)
      decompress   :  120.82 MVx/sec
      decompress+gz:  117.74 MVx/sec
    
RANDOM NOISE [0,2000) uint32

      compress     :  48.85 MVx/sec (8411876 bytes, 50.1%)
      compress+gz  :  25.54 MVx/sec (6815937 bytes, 81.0%)
      decompress   :  132.24 MVx/sec
      decompress+gz:  88.76 MVx/sec
    

      compress     :  48.22 MVx/sec (8411863 bytes, 50.1%)
      compress+gz  :  25.36 MVx/sec (6816509 bytes, 81.0%)
      decompress   :  135.53 MVx/sec
      decompress+gz:  90.15 MVx/sec
    

      compress     :  51.76 MVx/sec (8411476 bytes, 50.1%)
      compress+gz  :  25.67 MVx/sec (6816138 bytes, 81.0%)
      decompress   :  131.40 MVx/sec
      decompress+gz:  88.67 MVx/sec
    
BINARY NOISE [0,1] uint8 (pathological case)

      compress     :  26.93 MVx/sec (1820793 bytes, 43.4%)
      compress+gz  :  23.68 MVx/sec (1143454 bytes, 62.8%)
      decompress   :  61.18 MVx/sec
      decompress+gz:  58.87 MVx/sec
    

      compress     :  26.93 MVx/sec (1821063 bytes, 43.4%)
      compress+gz  :  23.73 MVx/sec (1143784 bytes, 62.8%)
      decompress   :  58.38 MVx/sec
      decompress+gz:  56.24 MVx/sec
    

      compress     :  27.89 MVx/sec (1822417 bytes, 43.4%)
      compress+gz  :  24.52 MVx/sec (1144408 bytes, 62.8%)
      decompress   :  58.50 MVx/sec
      decompress+gz:  56.40 MVx/sec
    
EMPTY

      compress     :  510.20 MVx/sec (1262 bytes, 0.0%)
      compress+gz  :  508.78 MVx/sec (77 bytes, 6.1%)
      decompress   :  5848.47 MVx/sec
      decompress+gz:  5823.30 MVx/sec
    
SOLID ONES

      compress     :  488.35 MVx/sec (1262 bytes, 0.0%)
      compress+gz  :  487.66 MVx/sec (79 bytes, 6.3%)
      decompress   :  4079.82 MVx/sec
      decompress+gz:  4063.80 MVx/sec
