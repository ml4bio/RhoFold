#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
rootdir=`dirname $bindir`
cd $rootdir

echo "download RNAcentral"
wget -q ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_species_specific_ids.fasta.gz -O rnacentral_species_specific_ids.fasta.gz 
wget -q ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/release_notes.txt -O release_notes.txt
wget -q ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/rfam/rfam_annotations.tsv.gz -O rfam_annotations.tsv.gz

echo "download rfam"
wget -q ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.full_region.gz -O Rfam.full_region.gz 
wget -q ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz -O Rfam.cm.gz 
wget -q ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/version.txt.gz -O version.txt.gz
gzip -d -f version.txt.gz

echo "download nt"
wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz -O nt.gz
wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz.md5 -O nt.gz.md5
md5sum -c nt.gz.md5
wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/db/ -O nt.html
for url in `grep -ohP "ftp://ftp.ncbi.nlm.nih.gov[:\d]*/blast/db/nt.\d+.tar.gz" nt.html |uniq`;do
    filename=`basename $url`
    wget -q $url -O $filename
    wget -q $url.md5 -O $filename.md5
    md5sum -c $filename.md5
done
rm *md5
