#!/usr/bin/env zsh
git commit -m "${1:-release}" main.nf
git push origin main
# add one to last number in current version
next_vers=`gh release list | grep Latest | perl -ne '($a) = /(v\d+\.\d+\.\d+)/; $a =~ s/(\d+)$/1+$1/e; print "$a\n" '`
echo Making next release: $next_vers
gh release create $next_vers  --notes "${1:-release}"
./run_job_as_batches.sh $next_vers

