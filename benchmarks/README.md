Note that all numbers provided are for the "flat" label layout. The pin layout is still under investigation (it is much slower and may not always provide a benefit). All data were first converted to Fortran order.

## Size Benchmark
shape: (256, 256, 64)  
PINKY40 CUTOUTS (connectomics.npy)

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
BINARY IMAGE (connectomics.npy, label 67699431)

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
WATERSHED CUTOUTS (ws.npy)

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
EMPTY

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
SOLID ONES

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
RANDOM NOISE [0,2000) uint32

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
BINARY NOISE [0,1] uint8 (pathological case)

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
PINKY40 CUTOUTS (connectomics.npy)

      compress     :  200.48 MVx/sec (32155 bytes, 0.2%)
      compress+gz  :  194.84 MVx/sec (23329 bytes, 72.6%)
      decompress   :  122.70 MVx/sec
      decompress+gz:  122.44 MVx/sec


      compress     :  236.59 MVx/sec (10900 bytes, 0.1%)
      compress+gz  :  234.43 MVx/sec (4647 bytes, 42.6%)
      decompress   :  124.95 MVx/sec
      decompress+gz:  124.88 MVx/sec


      compress     :  240.57 MVx/sec (11332 bytes, 0.1%)
      compress+gz  :  238.25 MVx/sec (5345 bytes, 47.2%)
      decompress   :  124.34 MVx/sec
      decompress+gz:  124.26 MVx/sec

WATERSHED CUTOUTS (ws.npy)

      compress     :  54.78 MVx/sec (568444 bytes, 1.7%)
      compress+gz  :  50.40 MVx/sec (419719 bytes, 73.8%)
      decompress   :  77.88 MVx/sec
      decompress+gz:  76.26 MVx/sec


      compress     :  62.43 MVx/sec (535607 bytes, 1.6%)
      compress+gz  :  57.23 MVx/sec (397003 bytes, 74.1%)
      decompress   :  80.39 MVx/sec
      decompress+gz:  78.75 MVx/sec


      compress     :  63.48 MVx/sec (547982 bytes, 1.6%)
      compress+gz  :  58.05 MVx/sec (402553 bytes, 73.5%)
      decompress   :  79.58 MVx/sec
      decompress+gz:  77.94 MVx/sec

RANDOM NOISE [0,2000) uint32

      compress     :  21.66 MVx/sec (8411095 bytes, 50.1%)
      compress+gz  :  14.34 MVx/sec (6815999 bytes, 81.0%)
      decompress   :  93.41 MVx/sec
      decompress+gz:  64.41 MVx/sec


      compress     :  20.56 MVx/sec (8411303 bytes, 50.1%)
      compress+gz  :  13.93 MVx/sec (6815973 bytes, 81.0%)
      decompress   :  95.64 MVx/sec
      decompress+gz:  66.80 MVx/sec


      compress     :  21.60 MVx/sec (8411656 bytes, 50.1%)
      compress+gz  :  14.37 MVx/sec (6816264 bytes, 81.0%)
      decompress   :  91.83 MVx/sec
      decompress+gz:  64.74 MVx/sec

BINARY NOISE [0,1] uint8 (pathological case)

      compress     :  18.31 MVx/sec (1822410 bytes, 43.4%)
      compress+gz  :  16.37 MVx/sec (1141655 bytes, 62.6%)
      decompress   :  39.07 MVx/sec
      decompress+gz:  37.94 MVx/sec


      compress     :  18.13 MVx/sec (1821937 bytes, 43.4%)
      compress+gz  :  16.22 MVx/sec (1141636 bytes, 62.7%)
      decompress   :  38.66 MVx/sec
      decompress+gz:  37.42 MVx/sec


      compress     :  18.26 MVx/sec (1821654 bytes, 43.4%)
      compress+gz  :  16.33 MVx/sec (1141202 bytes, 62.6%)
      decompress   :  39.22 MVx/sec
      decompress+gz:  38.06 MVx/sec

EMPTY

      compress     :  231.11 MVx/sec (993 bytes, 0.0%)
      compress+gz  :  230.73 MVx/sec (62 bytes, 6.2%)
      decompress   :  121.75 MVx/sec
      decompress+gz:  121.74 MVx/sec

SOLID ONES

      compress     :  250.69 MVx/sec (993 bytes, 0.0%)
      compress+gz  :  250.12 MVx/sec (66 bytes, 6.6%)
      decompress   :  122.90 MVx/sec
      decompress+gz:  122.89 MVx/sec