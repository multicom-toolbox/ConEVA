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
distance :19@CB :47@CB
distance :63@CA :108@CB
distance :53@CB :91@CB
distance :53@CB :112@CB
distance :19@CB :44@CB
distance :18@CB :88@CA
distance :18@CB :47@CB
distance :18@CB :46@CB
distance :53@CB :87@CB
distance :51@CB :88@CA
distance :25@CB :98@CB
distance :33@CA :72@CB
distance :4@CB :134@CB
distance :21@CB :79@CB
distance :55@CA :95@CB
distance :5@CA :73@CB
distance :54@CB :119@CB
distance :9@CB :71@CB
distance :39@CB :100@CA
distance :9@CB :138@CB
distance :53@CB :89@CB
distance :22@CB :85@CB
distance :57@CB :129@CB
distance :24@CB :79@CB
distance :18@CB :52@CB
distance :74@CB :105@CB
distance :6@CB :134@CB
distance :56@CB :93@CB
distance :4@CB :76@CB
distance :57@CB :94@CB
distance :61@CB :85@CB
distance :76@CB :100@CA
distance :39@CB :75@CA
distance :53@CB :119@CB
distance :58@CB :85@CB
distance :58@CB :125@CB
distance :56@CB :95@CB
distance :22@CB :52@CB
distance :59@CB :98@CB
distance :52@CB :112@CB
distance :63@CA :87@CB
distance :71@CB :119@CB
distance :58@CB :97@CB
distance :58@CB :112@CB
distance :35@CB :102@CB
distance :64@CB :103@CB
distance :61@CB :109@CB
distance :53@CB :94@CB
distance :74@CB :101@CB
distance :76@CB :102@CB
distance :43@CB :128@CB
distance :8@CB :75@CA
distance :54@CB :116@CB
distance :11@CB :71@CB
distance :52@CB :85@CB
distance :54@CB :95@CB
distance :4@CB :78@CB
distance :61@CB :87@CB
distance :58@CB :99@CB
distance :95@CB :141@CB
distance :9@CB :77@CB
distance :10@CB :79@CB
distance :54@CB :109@CB
distance :54@CB :120@CB
distance :9@CB :78@CB
distance :11@CB :78@CB
distance :75@CA :100@CA
distance :59@CB :103@CB
distance :5@CA :32@CA
distance :53@CB :85@CB
distance :71@CB :115@CB
distance :71@CB :103@CB
distance :111@CB :141@CB
distance :9@CB :76@CB
distance :30@CB :71@CB
distance :60@CB :108@CB
distance :71@CB :101@CB
distance :25@CB :52@CB
distance :35@CB :71@CB
distance :77@CB :109@CB
distance :64@CB :108@CB
distance :57@CB :125@CB
distance :10@CB :77@CB
distance :57@CB :89@CB
distance :54@CB :98@CB
distance :74@CB :102@CB
distance :78@CB :128@CB
distance :18@CB :51@CB
distance :55@CA :116@CB
distance :43@CB :115@CB
distance :35@CB :72@CB
distance :95@CB :128@CB
distance :13@CB :82@CB
distance :5@CA :34@CB
distance :13@CB :81@CB
distance :64@CB :116@CB
distance :58@CB :100@CA
distance :60@CB :87@CB
distance :78@CB :111@CB
distance :51@CB :87@CB
distance :51@CB :90@CB
distance :72@CB :104@CB
distance :13@CB :80@CB
distance :76@CB :103@CB
distance :115@CB :140@CB
distance :18@CB :48@CB
distance :77@CB :101@CB
distance :58@CB :87@CB
distance :57@CB :85@CB
distance :22@CB :61@CB
distance :71@CB :99@CB
distance :4@CB :111@CB
distance :60@CB :85@CB
distance :51@CB :84@CA
distance :26@CA :87@CB
distance :80@CB :131@CB
distance :78@CB :102@CB
distance :76@CB :104@CB
distance :29@CB :76@CB
distance :22@CB :103@CB
distance :62@CB :108@CB
distance :66@CB :105@CB
distance :53@CB :95@CB
distance :34@CB :72@CB
distance :25@CB :60@CB
distance :98@CB :128@CB
distance :55@CA :102@CB
distance :42@CB :69@CB
distance :43@CB :79@CB
distance :74@CB :111@CB
distance :29@CB :115@CB
distance :61@CB :100@CA
distance :75@CA :103@CB
distance :100@CA :128@CB
distance :62@CB :101@CB
distance :60@CB :113@CB
distance :24@CB :81@CB
distance :52@CB :103@CB
distance :50@CB :112@CB
distance :76@CB :120@CB
distance :60@CB :100@CA
distance :4@CB :138@CB
distance :74@CB :120@CB
distance :8@CB :77@CB
distance :11@CB :140@CB
distance :53@CB :120@CB
distance :61@CB :108@CB
distance :55@CA :103@CB
distance :5@CA :76@CB
distance :50@CB :101@CB
distance :16@CB :81@CB
distance :22@CB :74@CB
distance :74@CB :108@CB
distance :35@CB :70@CB
distance :74@CB :100@CA
distance :55@CA :89@CB
distance :76@CB :113@CB
distance :8@CB :98@CB
distance :75@CA :101@CB
distance :52@CB :94@CB
distance :22@CB :59@CB
distance :62@CB :87@CB
distance :34@CB :74@CB
distance :57@CB :104@CB
distance :8@CB :141@CB
distance :36@CB :73@CB
distance :11@CB :115@CB
distance :57@CB :93@CB
distance :22@CB :77@CB
distance :54@CB :112@CB
distance :25@CB :54@CB
distance :42@CB :71@CB
distance :76@CB :119@CB
distance :61@CB :99@CB
distance :65@CB :120@CB
distance :4@CB :97@CB
distance :37@CB :72@CB
distance :18@CB :42@CB
distance :43@CB :102@CB
distance :5@CA :132@CB
distance :55@CA :91@CB
distance :44@CB :78@CB
distance :58@CB :89@CB
distance :92@CB :118@CB
distance :25@CB :97@CB
distance :5@CA :138@CB
distance :79@CB :108@CB
distance :43@CB :103@CB
distance :85@CB :111@CB
distance :58@CB :101@CB
distance :33@CA :75@CA
distance :60@CB :120@CB
distance :11@CB :111@CB
distance :77@CB :112@CB
distance :29@CB :112@CB
distance :43@CB :77@CB
distance :12@CA :132@CB
distance :60@CB :102@CB
distance :54@CB :118@CB
distance :54@CB :108@CB
distance :34@CB :102@CB
distance :9@CB :54@CB
distance :4@CB :119@CB
distance :59@CB :99@CB
distance :54@CB :90@CB
distance :78@CB :131@CB
distance :78@CB :115@CB
distance :97@CB :140@CB
distance :53@CB :101@CB
distance :103@CB :141@CB
distance :17@CB :132@CB
distance :10@CB :134@CB
distance :52@CB :80@CB
distance :65@CB :103@CB
distance :51@CB :101@CB
distance :18@CB :85@CB
distance :29@CB :78@CB
setattr g lineWidth 2;
setattr g color black;
turn z 2 360 precess 20;
