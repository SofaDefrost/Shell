awk 'BEGIN { print "<Mesh position=\""; section=0 } 
/Coordinates/ { section=section+1 } 
/Triangle/ {section=2}
(NF==4&&section==1) { print $2 " " $3 " " $4 } 
/Elements/ {  
	     if (section==1) { print "\"\n tetrahedra=\""  } 
	     if (section==2)  {print "\"\n triangles=\"" } 
	    } 

(NF==5&&section==1) { print $2-1 " " $3-1 " " $4-1 " " $5-1 } 
(NF==4&&section==2) {print $2-1 " " $4-1 " " $3-1 }
END { print "\" />" }'
