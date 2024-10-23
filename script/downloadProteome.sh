#!/bin/bash

FASTADIR=/home/seamustard52/bfvd-analysis/proteome_cover/fasta

wget -O $FASTADIR/sars_cov_2.fasta "https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=%28%28proteome%3AUP000464024%29%29"
wget -O $FASTADIR/bacteriophage_t4.fasta "https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=%28%28proteome%3AUP000009087%29%29"
wget -O $FASTADIR/hhv_1.fasta "https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=%28%28proteome%3AUP000009294%29+AND+reviewed%3Dtrue%29"
wget -O $FASTADIR/tmv_u1.fasta "https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=%28%28proteome%3AUP000000522%29+AND+reviewed%3Dtrue%29"
wget -O $FASTADIR/sirv_1.fasta "https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=%28%28proteome%3AUP000002270%29%29"
wget -O $FASTADIR/hiv_1.fasta "https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=%28%28proteome%3AUP000007420%29+AND+reviewed%3Dtrue%29"
wget -O $FASTADIR/influenza_a.fasta "https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=%28%28proteome%3AUP000131152%29%29"
