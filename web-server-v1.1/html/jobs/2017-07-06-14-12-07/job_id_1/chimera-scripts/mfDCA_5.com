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
distance :43@CB :100@CA
distance :25@CB :100@CA
distance :75@CA :102@CB
distance :39@CB :102@CB
setattr g lineWidth 2;
setattr g color black;
turn z 2 360 precess 20;
