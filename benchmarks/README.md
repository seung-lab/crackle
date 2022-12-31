Note that all numbers provided are for the "flat" label layout. The pin layout is still under investigation (it is much slower and may not always provide a benefit). All data were first converted to Fortran order.

## Size Benchmark
shape: (256, 256, 64)  
PINKY40 CUTOUTS (connectomics.npy)

      ckl:         28348    
      cpso:        67864  (41.77%) 
      raw:      16777216  (0.17%)
      ckl.gz:      19228   
      cpso.gz      26256  (73.23%)
      raw:.gz      71158  (27.02%)
    

      ckl:         14516    
      cpso:        28496  (50.94%) 
      raw:      16777216  (0.09%)
      ckl.gz:       7727   
      cpso.gz       8963  (86.21%)
      raw:.gz      28636  (26.98%)
    

      ckl:         24740    
      cpso:        56894  (43.48%) 
      raw:      16777216  (0.15%)
      ckl.gz:      15809   
      cpso.gz      21184  (74.63%)
      raw:.gz      59809  (26.43%)
    
WATERSHED CUTOUTS (ws.npy)

      ckl:        537617    
      cpso:      2515670  (21.37%) 
      raw:      33554432  (1.60%)
      ckl.gz:     400151   
      cpso.gz     642979  (62.23%)
      raw:.gz     970262  (41.24%)
    

      ckl:        561911    
      cpso:      2678712  (20.98%) 
      raw:      33554432  (1.67%)
      ckl.gz:     421943   
      cpso.gz     680361  (62.02%)
      raw:.gz    1012515  (41.67%)
    

      ckl:        568996    
      cpso:      2670264  (21.31%) 
      raw:      33554432  (1.70%)
      ckl.gz:     428341   
      cpso.gz     694560  (61.67%)
      raw:.gz    1026001  (41.75%)
    
RANDOM NOISE [0,2000) uint32

      ckl:       8405450    
      cpso:     17301994  (48.58%) 
      raw:      16777216  (50.10%)
      ckl.gz:    6791314   
      cpso.gz    8024443  (84.63%)
      raw:.gz    8023562  (84.64%)
    

      ckl:       8405384    
      cpso:     17301996  (48.58%) 
      raw:      16777216  (50.10%)
      ckl.gz:    6791995   
      cpso.gz    8024532  (84.64%)
      raw:.gz    8023995  (84.65%)
    

      ckl:       8405714    
      cpso:     17301984  (48.58%) 
      raw:      16777216  (50.10%)
      ckl.gz:    6791443   
      cpso.gz    8024041  (84.64%)
      raw:.gz    8023287  (84.65%)
    
BINARY NOISE [0,1] uint8 (pathological case)

      ckl:       1842719    
      cpso:      2899167  (63.56%) 
      raw:       4194304  (43.93%)
      ckl.gz:    1193828   
      cpso.gz    1098040  (108.72%)
      raw:.gz     667063  (178.97%)
    

      ckl:       1843897    
      cpso:      2900235  (63.58%) 
      raw:       4194304  (43.96%)
      ckl.gz:    1194222   
      cpso.gz    1098982  (108.67%)
      raw:.gz     667228  (178.98%)
    

      ckl:       1843267    
      cpso:      2897477  (63.62%) 
      raw:       4194304  (43.95%)
      ckl.gz:    1195040   
      cpso.gz    1097542  (108.88%)
      raw:.gz     667036  (179.16%)
    
EMPTY

      ckl:           608    
      cpso:          824  (73.79%) 
      raw:      16777216  (0.00%)
      ckl.gz:         35   
      cpso.gz         43  (81.40%)
      raw:.gz      16316  (0.21%)
    
SOLID ONES

      ckl:           608    
      cpso:          824  (73.79%) 
      raw:      16777216  (0.00%)
      ckl.gz:         37   
      cpso.gz         47  (78.72%)
      raw:.gz      16322  (0.23%)
    


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
    
