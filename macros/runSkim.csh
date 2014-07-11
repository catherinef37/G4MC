#!/bin/csh -fb 

if ($#argv < 2 ) then
echo "usage: $0:t <level(1-8)> <rootfile_list>"
echo example: $0:t 1 nt1.root nt2.root nt3.root
exit -1
endif

set level = ($1)
foreach file ($argv[2-$#argv])
    set filename = ($file:r)
    ln -sf $file histfile.root
    root -b -l -q ../../macros/Skim.C\($level,\"track0\"\)
    mv histfile_skimmed.root ${filename}_skimmed${level}.root
    rm -f histfile.root
    ls -lh $file ${filename}_skimmed${level}.root
end
