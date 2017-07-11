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
distance :75@CA :102@CB
distance :71@CB :102@CB
distance :25@CB :79@CB
distance :21@CB :98@CB
distance :96@CB :129@CB
distance :76@CB :111@CB
distance :25@CB :77@CB
distance :43@CB :100@CA
distance :74@CB :103@CB
distance :12@CA :80@CB
distance :25@CB :100@CA
distance :39@CB :102@CB
distance :78@CB :137@CB
distance :29@CB :77@CB
distance :13@CB :79@CB
distance :34@CB :75@CA
distance :76@CB :101@CB
distance :9@CB :134@CB
distance :14@CB :82@CB
distance :80@CB :137@CB
distance :11@CB :80@CB
distance :9@CB :137@CB
distance :21@CB :85@CB
distance :29@CB :75@CA
distance :74@CB :104@CB
distance :38@CB :70@CB
distance :58@CB :119@CB
distance :15@CB :81@CB
distance :21@CB :81@CB
setattr g lineWidth 2;
setattr g color black;
turn z 2 360 precess 20;
