# images with blast, ngmlr, bcftools, samtools, bin_reads_by_umi.py

# shorthand for docker image pull
FROM ubuntu:22.04 

# first update, then install pretty generic tools
#==================================================
# python
# aws cli (v2)
# ncbi-blast+
# include prereqs for samtools
RUN apt-get update && apt -y install make cmake ncbi-blast+ curl bzip2 unzip gawk gcc g++ python3 python3-pip  tar wget libncurses-dev libbz2-dev liblzma-dev zlib1g-dev git

# removed bioperl
# had it in list for bp_gff2genbank.pl # this is way too big ~ 1Gb! do something else

# python libraries
#===================
RUN pip install numpy pandas matplotlib biopython

# AWS CLI
#==========
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install




# bioinf stuff
#===============

# my own python script and libs
COPY bin_reads_by_umi.py gb_parsing.py util.py Fastaq.py scatter.py histogram.py /usr/local/bin/

# samtools, installs into /usr/local/bin/
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && tar -xjf samtools-1.18.tar.bz2  && cd samtools-1.18 && ./configure && make && make install 
# use bcftools csq instead of snpeff
RUN wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 && tar -xjf bcftools-1.18.tar.bz2  && cd bcftools-1.18 && ./configure && make && make install
RUN wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 && tar -xjf htslib-1.18.tar.bz2 && cd htslib-1.18 && ./configure && make && make install

# ngmlr
# the live release 0.2.7 has a signed integer overflow error
# so clone from git and compile from source
# use this commit: a2a31fb which has this fix in src/
# Fix for "Signed Integer Overflow in MAPQ"
RUN git clone https://github.com/philres/ngmlr.git && cd ngmlr && mkdir build && cd build && cmake .. && make && cd ../bin && ln -s `pwd`/ngmlr-0.2.8/ngmlr  /usr/local/bin/ngmlr

# not using these pieces
#==========================
#CMD ["ngmlr"]
#ENTRYPOINT ["/bin/bash" , "-c"]  # breaks the image when run in NF
