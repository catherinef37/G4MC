#!/bin/csh -fb
#try to calculate how many events hit the VD
#usage $0:t <filelist>
 
foreach file ($*)
    if !( -f $file ) continue
    root -b $file <<EOF
    track0->Draw("Pvb>>h","Pvb>0");
    printf("file=%40s,  NHits=%5d\n",_file0->GetName(),h->GetEntries());
    .q
EOF

end
