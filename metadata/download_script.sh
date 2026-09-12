#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=24:00:00
#SBATCH --account=project_2005718
#SBATCH --mail-type=ALL
#SBATCH --partition=small
# Download LIFEPLAN new sequencing data from Dropbox
set -e
mkdir -p LIFEPLAN-00042 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/AO4xIsiU0qsEWFy93hxYq2U/LPLAN42_NSEQ-00103/CBGMB01554-1565_BF3-BR2_NSEQ-00103_S1_L001_R1_001.fastq.gz?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00042/CBGMB01554-1565_BF3-BR2_NSEQ-00103_S1_L001_R1_001.fastq.gz
mkdir -p LIFEPLAN-00042 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/AEa-nUqd9J4yR7dtr6I0B94/LPLAN42_NSEQ-00103/CBGMB01554-1565_BF3-BR2_NSEQ-00103_S1_L001_R2_001.fastq.gz?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00042/CBGMB01554-1565_BF3-BR2_NSEQ-00103_S1_L001_R2_001.fastq.gz
mkdir -p LIFEPLAN-00042 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/AKjnNr4SV4YuUw_JpQnfi6g/LPLAN42_NSEQ-00103/NSEQ-00103_ClarityLimsExtract.xlsx?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00042/NSEQ-00103_ClarityLimsExtract.xlsx
mkdir -p LIFEPLAN-00042 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/AO3L-DdIHgIpACZ6NdjarOM/LPLAN42_NSEQ-00103/NSEQ-00103_UMIMap.txt?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00042/NSEQ-00103_UMIMap.txt
mkdir -p LIFEPLAN-00043 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/AMFJ4Pw2vlO58vOiNL85jO8/LPLAN43_NSEQ-00090/CBGMB01396-1407_BF3-BR2_NSEQ-00090_S1_L001_R1_001.fastq.gz?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00043/CBGMB01396-1407_BF3-BR2_NSEQ-00090_S1_L001_R1_001.fastq.gz
mkdir -p LIFEPLAN-00043 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/ACa-7vh01TYdixcI6Mz9RyE/LPLAN43_NSEQ-00090/CBGMB01396-1407_BF3-BR2_NSEQ-00090_S1_L001_R2_001.fastq.gz?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00043/CBGMB01396-1407_BF3-BR2_NSEQ-00090_S1_L001_R2_001.fastq.gz
mkdir -p LIFEPLAN-00043 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/AJDA5r6t0s_GlwvNyjvSV_k/LPLAN43_NSEQ-00090/NSEQ-00090_ClarityLimsExtract.xlsx?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00043/NSEQ-00090_ClarityLimsExtract.xlsx
mkdir -p LIFEPLAN-00043 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/ANTB2OxYIs4oyUpzeAwwoC4/LPLAN43_NSEQ-00090/NSEQ-00090_UMIMap.txt?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00043/NSEQ-00090_UMIMap.txt
mkdir -p LIFEPLAN-00044 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/AGukDDzbh0d9WHne-ieGivs/LPLAN44_NSEQ-00094/CBGMB01444-1455_BF3-BR2_NSEQ-00094_S1_L001_R1_001.fastq.gz?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00044/CBGMB01444-1455_BF3-BR2_NSEQ-00094_S1_L001_R1_001.fastq.gz
mkdir -p LIFEPLAN-00044 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/AEB0wxtmTdkBsv1dOV0ppCk/LPLAN44_NSEQ-00094/CBGMB01444-1455_BF3-BR2_NSEQ-00094_S1_L001_R2_001.fastq.gz?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00044/CBGMB01444-1455_BF3-BR2_NSEQ-00094_S1_L001_R2_001.fastq.gz
mkdir -p LIFEPLAN-00044 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/AHZC-rtMygGEwL7MEXDLay0/LPLAN44_NSEQ-00094/NSEQ-00094_ClarityLimsExtract.xlsx?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00044/NSEQ-00094_ClarityLimsExtract.xlsx
mkdir -p LIFEPLAN-00044 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/AIiMeqkQMonCfwiEC07fwJg/LPLAN44_NSEQ-00094/NSEQ-00094_UMIMap.txt?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00044/NSEQ-00094_UMIMap.txt
mkdir -p LIFEPLAN-00045 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/AIvrXo8IKLhUldxIELGhTE8/LPLAN45_NSEQ-00105/CBGMB01580-1591_BF3-BR2_NSEQ-00105_S1_L001_R1_001.fastq.gz?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00045/CBGMB01580-1591_BF3-BR2_NSEQ-00105_S1_L001_R1_001.fastq.gz
mkdir -p LIFEPLAN-00045 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/API-sv6R4D2zE33Iz8Di2nU/LPLAN45_NSEQ-00105/CBGMB01580-1591_BF3-BR2_NSEQ-00105_S1_L001_R2_001.fastq.gz?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00045/CBGMB01580-1591_BF3-BR2_NSEQ-00105_S1_L001_R2_001.fastq.gz
mkdir -p LIFEPLAN-00045 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/ANAuZQNStj7sMr5avb4jHBs/LPLAN45_NSEQ-00105/NSEQ-00105_ClarityLimsExtract.xlsx?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00045/NSEQ-00105_ClarityLimsExtract.xlsx
mkdir -p LIFEPLAN-00045 && wget "https://www.dropbox.com/scl/fo/8dfgtgzpb4170wc8j922u/ACm0dpCtqoToFOkwavEebqU/LPLAN45_NSEQ-00105/NSEQ-00105_UMIMap.txt?rlkey=ya89mlki2sgm1o2azs14vnw8p&dl=0" -O LIFEPLAN-00045/NSEQ-00105_UMIMap.txt
