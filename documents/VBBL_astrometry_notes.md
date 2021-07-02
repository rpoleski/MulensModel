To get the image positions from VBBL, following C++ code has to be run:

```c++
_sols *Images;
Mag=VBBL.BinaryMag(s, q, y1, y2, rho, accuracy, &Images);
for(_curve *c=Images->first; c; c=c->next) {
    for(_point *p=c->first; p; p=p->next) {
        // Here p->x1 and p->x2 give next point
    }
}
delete Images;
```

This code was written for version 1.X. Check if it is the same (maybe simpler) in 2.X.

