ssh -X sunam235@caucluster.rz.uni-kiel.de
1. basic linux
2. linux commands
3. Bioinformatics basic understanding
   
   - copy from one folder to another:
  
  block of code:
  ```sh
  cp source destination
  ```

  inline code:
  this is the command `cp`
  gunc plot -d /PATH/TO/YOUR/diamond_output/METABAT__#-contigs.diamond.progenomes_2.1.out -g /PATH/TO/YOUR/genes_calls/gene_counts.json
4. Bash script
For executing a process
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err
#SBATCH --partition=base
#SBATCH --reservation=biol217


5. Activates the conda environment
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

 cd /work_beegfs/sunam235/Metagenomics/assembly (path of directory)

 >>txt.file    


6. Coping to desktop
scp sunam235@caucluster.rz.uni-kiel.de:/work_beegfs/sunam235/Metagenomics/3_coassembly/*.fastg
coping to the desktop from a local terminal

So... from raw reads......

METAGENOMICS (pipeline and commands)
1. Quality control
1.1 fastqc
   Checks quality
   Input- raw reads (fastq files)
   output-fastqc report
 In a aloop
 ```sh
 cp for i in *.gz; do fastqc $i -o output_folder/; done
```
1.2 fastp
Trims and checks quality as well.
input-raw reads fastq files
output-fastp report (clean reads)
```sh
fastp -i ? -I ? -R ? -o ? -O ? -t 6 -q 20
```
2. Assembly (megahit tool)
   input-cleaned reads (fastq file format)
   output-contigs (fasta file format)-consensus region from aligned clean reads.
```sh
megahit -1 sample1_R1_clean.fastq.gz -1 sample2_R1_clean.fastq.gz -1 sample3_R1_clean.fastq.gz -2 sample1_R2_clean.fastq.gz -2 sample2_R2_clean.fastq.gz -2 sample3_R2_clean.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 -o /PATH/TO/3_coassembly/ -t 12
```
2.1 Contigs visualisation
contigs in fasta file format converted to fastg file format
```sh
megahit_toolkit contig2fastg 99 final.contigs.fa > final.contigs.fastg
 ```                 
visualized in Bandage (installed on the desktop)
2.2 Counting contigs (grep tool)
```sh
grep -c ">" final.contigs.fa
 ```
2.3 Assessment of quality of assemblies (metaquast tool)
checks quality of assembled contigs
```sh
metaquast -t 6 -o /PATH/TO/3_metaquast -m 1000 final.contigs.fa
```
3. Mapping (Bowtie2)
identifying which reads contribute to which contig-genomic coverage
The contigs are first formatted into a format acceptable by ANVI'O 
```sh
anvi-script-reformat-fasta ? -o ? --min-len 1000 --simplify-names --report-file name_conversion.txt
 ```
index the mapping reference fasta file
```sh
module load bowtie2
bowtie2-build contigs.anvio.fa contigs.anvio.fa.index
```
then align clean reads with contigs
```sh
module load bowtie2
bowtie2 --very-fast -x contigs.anvio.fa.index -1 /PATH/TO/sample1_R1_clean.fastq.gz -2 /PATH/TO/sample1_R2_clean.fastq.gz -S SAMPLE.sam
```
output-contigs in sam file format
4. Contigs database generation (anvi-gen tool)
first convert the contigs (sam file format) to bam file format
```sh
module load samtools
samtools view -bS ? > bam_file.bam
```
then generate the database
```sh
anvi-gen-contigs-database -f contigs.anvio.fa -o contigs.db -n 'biol217'
```
information entailed in the contigs database generated;
Compute k-mer frequencies for each contig
Soft-split contigs longer than 20,000 bp into smaller ones
Open reading frames using Prodigal, the bacterial and archaeal gene-finding program

5. HMM search on contigs
Hidden Markov Model
searches for specific genes with a known function in a larger dataset
it identifies hits (bacterial single-copy core gene collections) among your genes to the collections.
```sh
anvi-run-hmms -c contigs.db
```

6. Visualisation of contigs database
Done in the terminal to give us access to the NODES to visualise the contigs.db (with which we performed  hmm search) 
```sh
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-display-contigs-stats contigs.db
```

7. Binning with ANVI'O
 7.1 Preparation of contig data for binning.
  7.1.1 sorting and indexing bam files
   ```sh
   for i in *.bam; do anvi-init-bam $i -o "$i".sorted.bam; done
   ```
  7.1.2 creating  anvio profiles
   this profile stores sample-specific information about the contigs
   therefore gives  properties for each contig in a single sample, based on mapping results.
    ```sh
    anvi-profile -i YOUR_SORTED.bam -c contigs.db --output-dir OUTPUT_DIR
    ```
  7.1.3 ANVI'O Profile merging.
   the profiles are merged into one (that is; the 3 profiles you have for your 3 samples must be merged into 1 whole Anvio profile)
   this is done by overlapping all the profiles alongside the contigs database
   ```sh
   anvi-merge /PATH/TO/SAMPLE1/PROFILE.db /PATH/TO/SAMPLE2/PROFILE.db /PATH/TO/SAMPLE3/PROFILE.db -o /PATH/TO/merged_profiles -c /PATH/TO/contigs.db --enforce-hierarchical-clustering
   ```
7.2 Binning
   clustering of contigs together to form MAGs-Metagenome Assembled Genomes
   different tools can be used: Metabat2, binsanity, MaxBin2 among others.
  7.2.1 Binning with MetaBat2
    ```sh
   anvi-cluster-contigs -p /PATH/TO/merged_profiles/PROFILE.db -c /PATH/TO/contigs.db -C METABAT --driver metabat2 --just-do-it --log-file log-metabat2
    anvi-summarize -p /PATH/TO/merged_profiles/PROFILE.db -c /PATH/TO/contigs.db -o SUMMARY_METABAT -C METABAT
      ```
  7.2.2 Binning with MaxBin2
   ```sh
   anvi-cluster-contigs -p /PATH/TO/merged_profiles/PROFILE.db -c /PATH/TO/contigs.db -C MAXBIN2 --driver maxbin2 --just-do-it --log-file log-maxbin2
   anvi-summarize -p /PATH/TO/merged_profiles/PROFILE.db -c /PATH/TO/contigs.db -o SUMMARY_MAXBIN2 -C MAXBIN2
   ```
7.3 MAG quality estimation & visualization
   estimates genome completeness and contamination levels.
   completeness: high quality-above 90% and medium quality 50-90%
   used bins from the MetaBat2 binning strategy
   estimation
   ```sh
   anvi-estimate-genome-completeness -c /PATH/TO/contigs.db -p /PATH/TO/merged_profiles/PROFILE.db -C METABAT
   ```
   visualization
   ```sh
   module load gcc12-env/12.1.0
   module load miniconda3/4.12.0
   conda activate anvio-8
   anvi-interactive -p /PATH/TO/merged_profiles/PROFILE.db -c /PATH/TO/contigs.db -C YOUR_COLLECTION
  ```
7.4 Bin refinement
   can be done using;
   multiple binners
   reassembly (assemble the clean reads and the MAGs)
   chimera (different contigs binned together) detection with GUNC tool.
   first, create a folder called SUMMARY containing all summary files of the archaea bins created.
    ```sh
   module load gcc12-env/12.1.0
   module load miniconda3/4.12.0
   conda activate anvio-8
   anvi-summarize -p ? -c ? --list-collections
   anvi-summarize -c ? -p ? -C ? -o ? --just-do-it
    ```
   7.4.1 Chimera detection
   activate the gunc environment
     ```sh
   module load gcc12-env/12.1.0
   module load miniconda3/4.12.0
   conda activate gunc
     ```
     then run the command
      ```sh
      cd /PATH/TO/ARCHAEA_BIN_REFINEMENT 
      mkdir GUNC
      for i in *.fa; do gunc run -i "$i" -r /work_beegfs/sunam###/Databases/gunc_db_progenomes2.1.dmnd --out_dir GUNC/"$i" --threads 10 --detailed_output; done
        ```
      
   then.........
   Run a QC on the MAGs again
  7.4.2 Manual bin refinement
   As large metagenome assemblies can result in hundreds of bins, pre-select the better ones for manual refinement, e.g. > 70% completeness.
   Before you start, make a copy/backup of your unrefined bins to avoid them from being overwritten.
   Use anvi refine to work on your bins manually.
   ```sh
   anvi-refine -c /PATH/TO/contigs.db -C METABAT -p /PATH/TO/merged_profiles/PROFILE.db --bin-id Bin_METABAT__##
   ```
   You can now sort your bins by GC content, by coverage or both.
   For refinement it is easier to use the clustering based on only differential coverage, and then only based on sequence composition in search for outliers.
   The interface allows you to categorize contigs into separate bins (selection tool). Unhighlighted contigs are removed when the data is saved.
   You can also evaluate taxonomy and duplicate single copy core genes.
   You can also remove contigs.
8. Classification of MAGs
      based on Average Nucleotide Identity (ANI) and Alignment fraction once the MAGs are aligned with referece sequnces from a genome database
      ANI: 95% and above- same species with reference organism
           less then 95% only the genus name is taken
           99% and above genus, species and strain name are taken
      Alignment fraction of 70% is the minimum considered for ANI to be calculated
      ```sh
      anvi-estimate-scg-taxonomy -c /PATH/TO/contigs.db -p /PATH/TO/profile.db --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy
      anvi-estimate-scg-taxonomy -c /PATH/TO/contigs.db -p /PATH/TO/profile.db --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy > temp.txt
      anvi-summarize -p /PATH/TO/merged_profiles/PROFILE.db -c /PATH/TO/contigs.db --metagenome-mode -o /PATH/TO/SUMMARY_METABAT2 -C METABAT2
      ```
GENOMICS (pipeline and commands)
   1. Quality control
      1.1 Short reads (fast qc & fatsp tool)
         ```sh
         #!/bin/bash
         #SBATCH --nodes=1
         #SBATCH --cpus-per-task=32
         #SBATCH --mem=128G
         #SBATCH --time=5:00:00
         #SBATCH --job-name=01_fastqc
         #SBATCH --output=01_fastqc.out
         #SBATCH --error=01_fastqc.err
         #SBATCH --partition=base
         #SBATCH --reservation=biol217
         module load gcc12-env/12.1.0
         module load miniconda3/4.12.0
         module load micromamba/1.4.2
         micromamba activate 01_short_reads_qc
         ## 1.1 fastqc raw reads
         # mkdir -p $WORK/genomics/1_short_reads_qc/1_fastqc_raw
         # for i in *.gz; do fastqc $i -o $WORK/genomics/1_short_reads_qc/1_fastqc_raw -t 32; done
         ## 1.2 fastp
         mkdir -p $WORK/genomics/1_short_reads_qc/2_cleaned_reads
         fastp -i $WORK/genomics/0_raw_reads/short_reads/241155E_R1.fastq.gz \
         -I $WORK/genomics/0_raw_reads/short_reads/241155E_R2.fastq.gz \
         -R $WORK/genomics/1_short_reads_qc/2_cleaned_reads/fastp_report \
         -h $WORK/genomics/1_short_reads_qc/2_cleaned_reads/report.html \
         -o $WORK/genomics/1_short_reads_qc/2_cleaned_reads/241155E_R1_clean.fastq.gz \
         -O $WORK/genomics/1_short_reads_qc/2_cleaned_reads/241155E_R2_clean.fastq.gz -t 6 -q 25

         jobinfo
         ```
         Check the quality again with fastqc
       1.2 long reads
         Nanoplot tool - checks quality
          ```sh
          micromamba activate 02_long_reads_qc

          cd $WORK/genomics/0_raw_reads/long_reads/
          NanoPlot --fastq $file -o $output_dir -t 6 --maxlength 40000 --minlength 1000 --plots kde --format png --N50 --dpi 300 --store --raw --tsv_stats --info_in_report
          ```
          filtlong tool- checks quality and trims
         ```sh
         filtlong --min_length 1000 --keep_percent 90 $file1 | gzip > sample1_cleaned_filtlong.fastq.gz
         mv sample1_cleaned_filtlong.fastq.gz $output_dir
         ```
         Check quality again with Nanoplot
      2. Hybrid assembly (Unicycler tool)
            Unicycler tool can clean the reads too if not cleaned prior
            assembles the clean Short reads with clean long reads
            ```sh
            micromamba activate 03_unicycler
            unicycler -1 $short_read1 -2 $short_read2 -l $long_reads -o $output_dir -t 32
            ```
      3. Quality assessment of assembly
               checks completeness and contamination levels
               3.1 Quast
               ```sh
               micromamba activate 04_checkm_quast
               
               quast.py assembly.fasta --circos -L --conserved-genes-finding --rna-finding\
               --glimmer --use-all-alignments --report-all-metrics -o $output_dir -t 16
               ```
               3.2 CheckM
               ```sh
               micromamba activate 04_checkm_quast
               # Create the output directory if it does not exist
               mkdir -p $checkm_out
               # Run CheckM for this assembly
               checkm lineage_wf $inputdir $output_dir -x fasta --tab_table --file $checkm_out/checkm_results -r -t 24
               # Run CheckM QA for this assembly
               checkm tree_qa ./$checkm_out
               checkm qa ./$checkm_out/lineage.ms ./$checkm_out/ -o 1 > ./$checkm_out/Final_table_01.csv
               checkm qa ./c$checkm_out/lineage.ms ./$checkm_out/ -o 2 > ./$checkm_out/final_table_checkm.csv
               ```
           3.3 Visualisation of Assemblies
               with Bandage on your desktop (if already installed)
         4. Assembly annotation (Prokka tool)
                  ```sh
                  micromamba activate 06_prokka
                  # Run Prokka on the file
                  prokka $input/assembly.fasta --outdir $output_dir --kingdom Bacteria --addgenes --cpus 32
                  ```
         5. Classification of genomes (GTDB-TK database)
                     first copy the .fna files of the annotated genomes to the GTDBTK directory
                     then...
                     ```sh
                     micromamba activate 07_gtdbtk
                     #run gtdb
                     gtdbtk classify_wf --cpus 12 --genome_dir $input_fna_files --out_dir $output_dir --extension .fna
                     #reduce cpu and increase the ram in bash script in order to have best performance
                     ```
            6. Combining reports (Multiqc tool)
                  multiqc will create the output directory on its own, so dont create it before running it
                  Run MultiQC to combine all the QC reports at once at the end of the pipeline.
                  ```sh
                  micromamba activate 01_short_reads_qc
                  # run multiqc
                  multiqc $input_dir -o $output_dir
                  ```
      PANGENOMICS
       1. Comparing genomes with ANVI'O
               Involves obtaining  consensus from comparing genomes of different organisms fro example strains of the same species.
               dendogram is drawn to show the relationship
               Usually done at a genus (here you will deal with a bigger number of genomes due to the diffeent species) or species level
               species level is better
               relevance: analyses genome function and functional relationships ANI can be computred and visualised
                Compared the genomes of 8 Vibrio jasicida strains using Anvio.
               information a pan-genome entails:
               multiple related genomes (singletons, accessory and core genes) and MAGsnumber of genes, gene cluster and orthologues (cluster: a group of genes with similar
               sequences
               Gene frequency analysis; analyses the frequency of; core, accessory and singleton genes
           1.1. Activate the conda environmemnt
           1.2 Download the data
               that is the  sequences of the different strains
               ```sh
               curl -L https://ndownloader.figshare.com/files/28965090 -o V_jascida_genomes.tar.gz
               tar -zxvf V_jascida_genomes.tar.gz
               ls V_jascida_genomes
               ```
           1.3 Creat contigs database from fasta files
                ```sh
               cd $WORK/pangenomics_test/V_jascida_genomes/
                ls *fasta | awk 'BEGIN{FS="_"}{print $1}' > genomes.txt

                # remove all contigs <2500 not
                for g in `cat genomes.txt`
                do
                echo
                echo "Working on $g ..."
                echo

                anvi-script-reformat-fasta ${g}_scaffolds.fasta \
                               --min-len 2500 \
                               --simplify-names \
                               -o ${g}_scaffolds_2.5K.fasta
                done

                # generate contigs.db
                for g in `cat genomes.txt`
                do
                echo
                echo "Working on $g ..."
                echo
                anvi-gen-contigs-database -f ${g}_scaffolds_2.5K.fasta \
                              -o V_jascida_${g}.db \
                              --num-threads 4 \
                              -n V_jascida_${g}
                done
                # annotate contigs.db

                for g in *.do
                do
                anvi-run-hmms -c $g --num-threads 4
                anvi-run-ncbi-cogs -c $g --num-threads 4
                anvi-scan-trnas -c $g --num-threads 4
                anvi-run-scg-taxyonomy -c $g --num-threads 4
                done
                ```
                1.4 Visualisation of contigs
                ```sh
               module load gcc12-env/12.1.0
               module load miniconda3/4.12.0
               conda activate anvio-8

               anvi-display-contigs-stats /path/to.your/databases/*db
                ```
                1.5 Create an external genomes file
                ```sh
               anvi-script-gen-genomes-file --input-dir /path/to/input/dir \
                             -o external-genomes.txt
                 ```
                1.6 Estimate contamination
                done in the terminal
                ```sh
                cd V_jascida_genomes
                anvi-estimate-genome-completeness -e external-genomes.txt
                ```
                estimated completeness too for some of my practice samples below,
               +----------------------------------------------------------------------------+----------+--------------+----------------+----------------+--------------+----------------+
| genome name                                                                | domain   |   confidence |   % completion |   % redundancy |   num_splits |   total length |
+============================================================================+==========+==============+================+================+==============+================+
| my_genome_Acinetobacter_pittii_PHEA-2                                      | BACTERIA |            1 |            100 |              0 |          193 |        3862530 |
+----------------------------------------------------------------------------+----------+--------------+----------------+----------------+--------------+----------------+
| my_genome_Bacillus-smithii-strain_DSM_4216                                 | BACTERIA |            1 |            100 |              0 |          168 |        3368778 |
+----------------------------------------------------------------------------+----------+--------------+----------------+----------------+--------------+----------------+
| my_genome_Bacillus_subtilis_subsp._subtilis_str._168                       | BACTERIA |            1 |            100 |              0 |          209 |        4215606 |
+----------------------------------------------------------------------------+----------+--------------+----------------+----------------+--------------+----------------+
| my_genome_Mycobacterium_tuberculosis_H37Rv                                 | BACTERIA |            1 |          98.59 |           1.41 |          220 |        4411532 |
+----------------------------------------------------------------------------+----------+--------------+----------------+----------------+--------------+----------------+
| my_genome_Pseudomonas_aeruginosa_strain_NCTC10332                          | BACTERIA |            1 |          97.18 |           1.41 |          315 |        6316979 |
+----------------------------------------------------------------------------+----------+--------------+----------------+----------------+--------------+----------------+
| my_genome_Salmonella_enterica_subsp._enterica_serovar_Typhimurium_str._LT2 | BACTERIA |            1 |          98.59 |              0 |          242 |        4857450 |
+----------------------------------------------------------------------------+----------+--------------+----------------+----------------+--------------+----------------+
| my_genome_Staphylococcus_aureus_subsp._aureus_NCTC_8325                    | BACTERIA |            1 |            100 |           2.82 |          141 |        2821361 |
+----------------------------------------------------------------------------+----------+--------------+----------------+----------------+--------------+----------------+
| my_genome_Staphylococcus_epidermidis-strain-ATCC_14990                     | BACTERIA |            1 |            100 |           1.41 |          123 |        2466502 |
+----------------------------------------------------------------------------+----------+--------------+----------------+----------------+--------------+----------------+
                1.7 Visualise contigs for refinement
               ```sh
               anvi-profile -c V_jascida_52.db \
               --sample-name V_jascida_52 \
               --output-dir V_jascida_52 \
               --blank
               ```
            1.8 Split genome into good bins
                ```sh
              anvi-split -p V_jascida_52/PROFILE.db \
              -c V_jascida_52.db \
               -C default \
               -o V_jascida_52_SPLIT
               # Here are the files you created
               #V_jascida_52_SPLIT/V_jascida_52_CLEAN/CONTIGS.db

               sed 's/V_jascida_52.db/V_jascida_52_SPLIT\/V_jascida_52_CLEAN\/CONTIGS.db/g' external-genomes.txt > external-genomes-final.txt
                 
   1.9 Estimate completeness of spit versus unsplit
     ```sh
   anvi-estimate-genome-completeness -e external-genomes.txt
                anvi-estimate-genome-completeness -e external-genomes-final.txt
                ```
         1.10 Compute the pangenome
               ```sh
               anvi-gen-genomes-storage -e external-genomes-final.txt \
                         -o V_jascida-GENOMES.db
anvi-pan-genome -g V_jascida-GENOMES.db \
                --project-name V_jascida \
                --num-threads 4                         
               ```
         1.11 Display the pangenome

           ```sh
         srun --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=base /bin/bash

                module load gcc12-env/12.1.0
                module load miniconda3/4.12.0
                conda activate anvio-8_biol217

                anvi-display-pan -p V_jascida/V_jascida-PAN.db \
                 -g V_jascida-GENOMES.db
                ```
   TRANSCRIPTOMICS
         1. RNAseq
               a transcriptome profiling technique used to identify and determine expression levels of transcripts
               involves:
               RNA-Seq data pre-processing
               RNA-Seq Data Analysis
               Differential Gene Expression
               Data Visualization
               SRA dataset used was from Prasse et al. 2017
         1.1 Quality control
               quality of the SRR sequences were checked using fastqc tool
               quality was good, so no trimming was done.
         1.2 READemption analysis
               READemption tool was used for the RNAseq analysis and visualisation was done in IGB
      
      ```sh
               #!/bin/bash
               #SBATCH --nodes=1
               #SBATCH --cpus-per-task=32
               #SBATCH --mem=64G
               #SBATCH --time=0-04:00:00
               #SBATCH --job-name=rna_seq_methanosarcina
               #SBATCH --output=rna_seq_methanosarcina.out
               #SBATCH --error=rna_seq_methanosarcina.err
               #SBATCH --partition=base
               #SBATCH --reservation=biol217

               module load gcc12-env/12.1.0
               module load miniconda3/4.12.0
               conda activate reademption

               ## 1. create a directory for the analysis
               reademption create --project_path READemption_analysis \
               --species metanosarcina="Methanosarcina mazei Gö1"

               #2- copy the sequences and files in respective directories
               # download the sequences from the NCBI database or github folder named "genome_input"

               #3- Processing and aligning the reads
               reademption align --project_path READemption_analysis \
               --processes 32 --segemehl_accuracy 95 \
               --poly_a_clipping \
               --fastq --min_phred_score 25 \
               --progress

               #4- Coverage
               reademption coverage --project_path READemption_analysis \
               --processes 32
               #5- Performing gene wise quantification
               reademption gene_quanti --project_path READemption_analysis \
               --processes 32 --features CDS,tRNA,rRNA
               #6- Performing differential gene expression analysis

               ####NOTE:: Change the names according to your file names in the READemption_analysis/input/reads/ directory

               reademption deseq --project_path READemption_analysis \
               --libs mut_R1.fastq.gz,mut_R2.fastq.gz,wt_R1.fastq.gz,wt_R2.fastq.gz \
               --conditions mut,mut,wt,wt --replicates 1,2,1,2 \
               --libs_by_species metanosarcina=mut_R1,mut_R2,wt_R1,wt_R2

               #7- Create plots
               reademption viz_align --project_path READemption_analysis
               reademption viz_gene_quanti --project_path READemption_analysis
               reademption viz_deseq --project_path READemption_analysis

               # The whole command will take around 2 hours to run.
               conda deactivate
               module purge
               jobinfo
               ```
   2. Riboprofiling (Riboseq)
                  study of ribosome-protected mRNA footprints
               can be used to determine translation efficiency, that is; the ratio of Ribo-Seq levels to RNA-seq levels
                  helps identify newly found genes and distinguish between coding and non-coding sequences.
                  results shoiuld be compared with other confirmatory analyses such as RNAseq.
                  SRR sequences were riboprofiled using HRIBO tool.
                visualisation of expression was done in IGB
                  
         R & R studio
