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
PINKY40 CUTOUTS (connectomics.npy) uint32

      compress     :  70.37 MVx/sec (151425 bytes, 0.9%)
      compress+gz  :  66.79 MVx/sec (114134 bytes, 75.4%)
      decompress   :  64.64 MVx/sec
      decompress+gz:  64.01 MVx/sec
    

      compress     :  158.80 MVx/sec (31232 bytes, 0.2%)
      compress+gz  :  156.48 MVx/sec (20619 bytes, 66.0%)
      decompress   :  78.27 MVx/sec
      decompress+gz:  78.12 MVx/sec
    

      compress     :  100.12 MVx/sec (96081 bytes, 0.6%)
      compress+gz  :  96.99 MVx/sec (71380 bytes, 74.3%)
      decompress   :  72.28 MVx/sec
      decompress+gz:  71.88 MVx/sec
    
WATERSHED CUTOUTS (ws.npy) uint64

      compress     :  29.90 MVx/sec (553901 bytes, 1.7%)
      compress+gz  :  28.50 MVx/sec (415444 bytes, 75.0%)
      decompress   :  48.43 MVx/sec
      decompress+gz:  47.46 MVx/sec
    

      compress     :  28.82 MVx/sec (584220 bytes, 1.7%)
      compress+gz  :  27.48 MVx/sec (434849 bytes, 74.4%)
      decompress   :  47.77 MVx/sec
      decompress+gz:  46.79 MVx/sec
    

      compress     :  30.49 MVx/sec (553152 bytes, 1.6%)
      compress+gz  :  29.05 MVx/sec (412760 bytes, 74.6%)
      decompress   :  47.81 MVx/sec
      decompress+gz:  46.87 MVx/sec
    
RANDOM NOISE [0,2000) uint32

      compress     :  17.97 MVx/sec (8405758 bytes, 50.1%)
      compress+gz  :  12.68 MVx/sec (6883083 bytes, 81.9%)
      decompress   :  58.81 MVx/sec
      decompress+gz:  41.49 MVx/sec
    

      compress     :  18.80 MVx/sec (8406123 bytes, 50.1%)
      compress+gz  :  13.09 MVx/sec (6883367 bytes, 81.9%)
      decompress   :  69.79 MVx/sec
      decompress+gz:  46.76 MVx/sec
    

      compress     :  18.79 MVx/sec (8405633 bytes, 50.1%)
      compress+gz  :  13.08 MVx/sec (6883017 bytes, 81.9%)
      decompress   :  69.37 MVx/sec
      decompress+gz:  46.53 MVx/sec
    
BINARY NOISE [0,1] uint8 (pathological case)

      compress     :  8.10 MVx/sec (1842845 bytes, 43.9%)
      compress+gz  :  7.69 MVx/sec (1193410 bytes, 64.8%)
      decompress   :  21.66 MVx/sec
      decompress+gz:  21.12 MVx/sec
    

      compress     :  8.12 MVx/sec (1843440 bytes, 44.0%)
      compress+gz  :  7.70 MVx/sec (1194535 bytes, 64.8%)
      decompress   :  21.77 MVx/sec
      decompress+gz:  21.22 MVx/sec
    

      compress     :  8.09 MVx/sec (1844092 bytes, 44.0%)
      compress+gz  :  7.68 MVx/sec (1193848 bytes, 64.7%)
      decompress   :  21.34 MVx/sec
      decompress+gz:  20.81 MVx/sec
    
EMPTY

      compress     :  171.12 MVx/sec (608 bytes, 0.0%)
      compress+gz  :  170.95 MVx/sec (47 bytes, 7.7%)
      decompress   :  76.28 MVx/sec
      decompress+gz:  76.27 MVx/sec
    
SOLID ONES

      compress     :  213.62 MVx/sec (608 bytes, 0.0%)
      compress+gz  :  213.39 MVx/sec (49 bytes, 8.1%)
      decompress   :  81.04 MVx/sec
      decompress+gz:  81.03 MVx/sec
    
