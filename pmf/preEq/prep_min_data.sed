#!/bin/bash

# Append Atom Names as comments to Masses Section
sed '/^Masses/,/^Pair Coeff/{
   s/^1 .*/& \# ODw/
   s/^2 .*/& \# H/
   s/^3 .*/& \# M/
   s/^4 .*/& \# O3m/
   s/^5 .*/& \# O3e/
   s/^6 .*/& \# polNa/
   s/^7 .*/& \# polI/
   s/^8 .*/& \# polBr/
}' data.slabPlumed.tip4p.$1.min > temp.$1.min

# Delete Pair Coeffs Section
sed '/^Pair Coeffs/,/^Bond Coeffs/{
  /^Bond Coeffs/ !{
     d
  }
}' temp.$1.min > data.slabPlumed.tip4p.$1.min

rm temp.$1.min
