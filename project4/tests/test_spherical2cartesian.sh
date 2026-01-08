#! /bin/bash
# This script tests the Spherical2Cartesian code which calculates the
# Transformation matrix between spherical and cartesian GTOs.

echo "Overlap matrix obtained using the MacMurchie-Davidson program:"

./../build/Spherical2Cartesian << EOF
3
EOF


echo -e "\n\n"
echo -e "  Expected transformation matrix for l = 3 for comparison and verification of correctness\n"
echo -e "          m :      -3          -2          -1           0           1           2           3         

x^3 y^0 z^0 :     0.000000    0.000000    0.000000    0.000000   -0.612372    0.000000    0.790569
x^0 y^3 z^0 :    -0.790569    0.000000   -0.612372    0.000000    0.000000    0.000000    0.000000
x^0 y^0 z^3 :     0.000000    0.000000    0.000000    1.000000    0.000000    0.000000    0.000000
x^2 y^1 z^0 :     2.371708    0.000000   -0.612372    0.000000    0.000000    0.000000    0.000000
x^2 y^0 z^1 :     0.000000    0.000000    0.000000   -1.500000    0.000000    1.936492    0.000000
x^1 y^2 z^0 :     0.000000    0.000000    0.000000    0.000000   -0.612372    0.000000   -2.371708
x^0 y^2 z^1 :     0.000000    0.000000    0.000000   -1.500000    0.000000   -1.936492    0.000000
x^1 y^0 z^2 :     0.000000    0.000000    0.000000    0.000000    2.449490    0.000000    0.000000
x^0 y^1 z^2 :     0.000000    0.000000    2.449490    0.000000    0.000000    0.000000    0.000000
x^1 y^1 z^1 :     0.000000    3.872983    0.000000    0.000000    0.000000    0.000000    0.000000
"
