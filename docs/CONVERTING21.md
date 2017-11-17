# Converting from older code for GooFit 2.1

The primary change in GooFit 2.1 was the addition of Python bindings, but there are several changes to watch out for. The ModernizeGooFit script will try to help, but it can't understand the details of your code. The options for the script have been expanded, check with `-h`. To run:

```
./scripts/ModernizeGooFit.py -v 2.1 filename.cu
```


## Variables and Observables

The three argument `Variable("name", min, max)` independent variable constructor has been separated out into a new class, `Observable`, which is now required everywhere for independent variables. The API should help you transition, since it's a compiler error to use one in place of the other.

Variables and observables are now passed by copy almost everywhere; they internally have shared pointers that keep each one synchronised.

## DecayInfo

The DecayInfo structs are now split up into four structs:

* `DecayInfo3`
* `DecayInfo3t`
* `DecayInfo4`
* `DecayInfo4t`

These are passed around by copy; but warning; they do hold raw pointers (explicitly handled by the python bindings).

## Resonances and Lineshapes

Resonances and Lineshapes are now separate constructors in `Resonances::` and `Lineshapes::`. The old enums were removed.



