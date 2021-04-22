
# ------------------------------------------------------------------------
if( length(list.files('data/hg19', 'PQS.hg19.bed', full.names = T)) == 0) {
  system(paste("python data/fastaRegexFinder.py", 
               "-f /mnt/0E471D453D8EE463/genomeDir/UCSC/hg19/hg19.fa",
               "-r '([gG]{3,}\\w{1,7}){3,}[gG]{3,}'",
               "> data/hg19/canonical_PQS_hg19.bed"))
  
  system(paste("python data/fastaRegexFinder.py", 
               "-f /mnt/0E471D453D8EE463/genomeDir/UCSC/hg19/hg19.fa",
               "-r '([gG]{2}\\w{1,12}){3,}[gG]{2}'",
               "> data/hg19/two_tetrads_PQS_hg19.bed"))
  
  system(paste("python data/fastaRegexFinder.py", 
               "-f /mnt/0E471D453D8EE463/genomeDir/UCSC/hg19/hg19.fa",
               "-r '([gG]{3,}\\w{8,12}){3,}[gG]{3,}'",
               "> data/hg19/extended_canonical_PQS_hg19.bed"))
}

# transtrand G4
# AAAA == canonical_PQS.bed
system(paste("python data/fastaRegexFinder.py", 
             "-f /mnt/0E471D453D8EE463/genomeDir/UCSC/hg19/hg19.fa",
             "-r '([gG]{3,}\\w{1,7}){3,}[gG]{3,}'",
             "> data/hg19/AAAA_hg19.bed"))

# AAAB
system(paste("python data/fastaRegexFinder.py", 
             "-f /mnt/0E471D453D8EE463/genomeDir/UCSC/hg19/hg19.fa",
             "-r '([gG]{3,}\\w{1,7}){3,}[cC]{3,}'",
             "> data/hg19/AAAB_hg19.bed"))

# AABA
system(paste("python data/fastaRegexFinder.py", 
             "-f /mnt/0E471D453D8EE463/genomeDir/UCSC/hg19/hg19.fa",
             "-r '([gG]{3,}\\w{1,7}){2,}[cC]{3,}\\w{1,7}[gG]{3,}'",
             "> data/hg19/AABA_hg19.bed"))

# ABAB
system(paste("python data/fastaRegexFinder.py", 
             "-f /mnt/0E471D453D8EE463/genomeDir/UCSC/hg19/hg19.fa",
             "-r '[gG]{3,}\\w{1,7}[cC]{3,}\\w{1,7}[gG]{3,}\\w{1,7}[cC]{3,}'",
             "> data/hg19/ABAB_hg19.bed"))

# ABAA
system(paste("python data/fastaRegexFinder.py", 
             "-f /mnt/0E471D453D8EE463/genomeDir/UCSC/hg19/hg19.fa",
             "-r '[gG]{3,}\\w{1,7}[cC]{3,}\\w{1,7}[gG]{3,}\\w{1,7}[gG]{3,}'",
             "> data/hg19/ABAA_hg19.bed"))

# ABBB
system(paste("python data/fastaRegexFinder.py", 
             "-f /mnt/0E471D453D8EE463/genomeDir/UCSC/hg19/hg19.fa",
             "-r '[gG]{3,}\\w{1,7}([cC]{3,}\\w{1,7}){3,}'",
             "> data/hg19/ABBB_hg19.bed"))

# ABBA
system(paste("python data/fastaRegexFinder.py", 
             "-f /mnt/0E471D453D8EE463/genomeDir/UCSC/hg19/hg19.fa",
             "-r '[gG]{3,}\\w{1,7}[cC]{3,}\\w{1,7}[cC]{3,}\\w{1,7}[gG]{3,}'",
             "> data/hg19/ABBA_hg19.bed"))

# AABB
system(paste("python data/fastaRegexFinder.py", 
             "-f /mnt/0E471D453D8EE463/genomeDir/UCSC/hg19/hg19.fa",
             "-r '([gG]{3,}\\w{1,7}){2,}([cC]{3,}\\w{1,7}){2,}'",
             "> data/hg19/AABB_hg19.bed"))
