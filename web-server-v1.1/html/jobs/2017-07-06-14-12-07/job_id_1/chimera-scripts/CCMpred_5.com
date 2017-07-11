# How to visualize these contacts in UCSF Chimera?
# 1. Download and install UCSF Chimera 1.10 or later if you do not have it already.
# 2. Save this file, say, 'all_sr.com'.
# 3. Use 'File->Open' in UCSF Chimera to open 'all_sr.com'.

background solid white;
open http://iris.rnet.missouri.edu/coneva/jobs/2017-07-06-14-12-07/job_id_1/native.pdb;
focus;
rainbow;
labelopt resinfo %(number)s
represent wire;
show;
rlabel;
focus;
distance :96@CB :129@CB
distance :6@CB :135@CB
distance :25@CB :100@CA
distance :14@CB :82@CB
distance :12@CA :80@CB
setattr g lineWidth 2;
setattr g color black;
turn z 2 360 precess 20;
