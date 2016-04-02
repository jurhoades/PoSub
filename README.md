# PoSub
A command line tool for finding codon point mutations that cause a specific amino acid substitution.

#### Inputs:
- Starting codon (3 letter string)
- Desired amino acid substitution (single letter, stop codon = "X")
- Optional: codon position for point mutation

A basic query:
```
python posub.py TCG L
```
A query with position argument:
```
python posub.py GAC E -p 3
```

Use the `-v` argument for more text output.

Tested for Python 2.7 and Python 3.5
