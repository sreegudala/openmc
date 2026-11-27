# N5 transposed k

```
TRIGA_DIR /home/sg44769/TRIGA_proposal_work/ce_mg_triga, TRIGA_LIB_DIR /home/sg44769/dev/triga_lib
N5: transposed library k on the TRIGA model against ce_mg_triga's forward k (kinetics_lib.md section 9)
  CE k 1.00080 +- 0.00050
  vector: transposed 1.00104 +- 0.00015 (20 replicas; spread 57 pcm, mean within-run sigma 67 pcm), forward 1.00137 +- 0.00009: transposed - forward -33 +- 17 pcm (-1.9 sigma) -> consistent; against CE +24 pcm
  matrix: transposed 1.00152 +- 0.00017 (20 replicas; spread 75 pcm, mean within-run sigma 66 pcm), forward 1.00152 +- 0.00009: transposed - forward -0 +- 19 pcm (-0.0 sigma) -> consistent; against CE +72 pcm
```
