#! /usr/bin/awk -f
# Converts MSH file (Gmsh) to WaveFront OBJ
#
# http://www.geuz.org/gmsh/doc/texinfo/gmsh.html#MSH-ASCII-file-format
#
# Important stuff:
#
#     $Nodes
#     number-of-nodes
#     node-number x-coord y-coord z-coord
#     ...
#     $EndNodes
#
#     $Elements
#     number-of-elements
#     elm-number elm-type number-of-tags < tag > ... node-number-list
#     ...
#     $EndElements
#
# NOTE: we assume the vertex numbers are ordered and continuous (no renumbering is done)
#

BEGIN { section = 0; }
/\$Nodes/ { section = 1; }
/\$EndNodes/ { section = 0; }
/\$Elements/ { section = 2; }
/\$EndElements/ { section = 0; }

# node-number x-coord y-coord z-coord
NF == 4 && section == 1 { print "v", $2, $3, $4; }

# elm-number elm-type number-of-tags < tag > ... node-number-list
# elm-type: 1 edges, 2 triangles
# ignore any tags!
NF == 6 && $2 == 2 && $3 == 0 && section == 2 { print "f", $4, $5, $6; }
