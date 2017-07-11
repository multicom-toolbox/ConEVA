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
distance :59@CB :95@CB
distance :97@CB :128@CB
distance :26@CA :100@CA
distance :12@CA :82@CB
distance :53@CB :90@CB
distance :34@CB :73@CB
distance :39@CB :71@CB
distance :29@CB :102@CB
distance :57@CB :96@CB
distance :54@CB :103@CB
distance :42@CB :70@CB
distance :11@CB :137@CB
distance :59@CB :87@CB
distance :28@CB :77@CB
distance :36@CB :72@CB
distance :59@CB :85@CB
distance :76@CB :115@CB
distance :97@CB :125@CB
distance :18@CB :87@CB
distance :6@CB :135@CB
distance :56@CB :125@CB
distance :42@CB :67@CB
distance :8@CB :34@CB
distance :97@CB :129@CB
distance :39@CB :63@CA
distance :6@CB :138@CB
distance :8@CB :32@CA
distance :60@CB :101@CB
distance :36@CB :70@CB
distance :18@CB :44@CB
distance :15@CB :79@CB
distance :11@CB :134@CB
distance :28@CB :79@CB
distance :7@CB :32@CA
distance :60@CB :112@CB
distance :78@CB :141@CB
distance :52@CB :87@CB
distance :29@CB :71@CB
distance :5@CA :75@CA
distance :14@CB :83@CB
distance :63@CA :102@CB
distance :77@CB :102@CB
distance :35@CB :75@CA
distance :74@CB :107@CB
setattr g lineWidth 2;
setattr g color black;
turn z 2 360 precess 20;
