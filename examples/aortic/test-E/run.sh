#! /bin/bash
cd $(dirname "$0")

for f in cyl1-E1e{4..7} cyl2-E1e{3..7} cyl2-bez-E1e{3..7} # example06-E1e{4,5,6,7}
do
    echo "===  $f ==="
    runSofa -l libshells.so -g batch -n 300 $f.scn
    ./msh2obj.awk < $f.gmsh > $f.obj
done
