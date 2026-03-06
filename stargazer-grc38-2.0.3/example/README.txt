# Get the path to the Stargazer directory

x=/YOUR/STARGAZER/DIRECTORY/PATH

# Move to the example directory

cd $x/example

# Display all of the Stargazer tools

python $x/stargazer --help

##############################################################################
# Genotype analysis for WGS data (GeT-RM; N=70) from Lee et al. (2019)       #
##############################################################################

# Below uses VDR as the control gene for copy number (CN) analysis

python $x/stargazer \
  -o getrm-cyp2d6-vdr \
  -d wgs \
  -t cyp2d6 \
  -c vdr \
  --vcf getrm-cyp2d6-vdr.joint.filtered.vcf \
  --gdf getrm-cyp2d6-vdr.gdf

# It's recommended to check CN results using multiple control genes

# Finally, you can provide a custom region as control

python $x/stargazer \
  -o getrm-cyp2d6-custom \
  -d wgs \
  -t cyp2d6 \
  --control_type custom \
  --region chr22:42546883-42551883 \
  --vcf getrm-cyp2d6-vdr.joint.filtered.vcf \
  --gdf getrm-cyp2d6-vdr.gdf

##############################################################################
# Genotype analysis for TS data (PGRNseq; N=96) from Lee et al. (2018)       #
##############################################################################

# Unlike WGS data, TS data requires inter-sample normalization for CN analysis
# Below uses the population mean during inter-sample normalization

python $x/stargazer \
  -o hapmap-cyp2d6-vdr \
  -d ts \
  -t cyp2d6 \
  -c vdr \
  --vcf hapmap-cyp2d6-vdr.joint.filtered.vcf \
  --gdf hapmap-cyp2d6-vdr.gdf

# For CN analysis, you may indicate known reference samples without SV
# Below will use the mean of indicated samples instead of the population mean

python $x/stargazer \
  -o hapmap-cyp2d6-vdr-ref \
  -d ts \
  -t cyp2d6 \
  -c vdr \
  --sample_list 133419 133420 133421 133423 133425 \
  --vcf hapmap-cyp2d6-vdr.joint.filtered.vcf \
  --gdf hapmap-cyp2d6-vdr.gdf

##############################################################################
# Genotype analysis with a VCF file only (WGS data)                          #
##############################################################################

python $x/stargazer \
  -o getrm-cyp3a5-vcfonly \
  -d wgs \
  -t cyp3a5 \
  --vcf getrm-cyp3a5-vdr.joint.filtered.vcf

##############################################################################
# Genotype analysis with a VCF file only (SNP array data)                    #
##############################################################################

python $x/stargazer \
  -o rok-cyp3a5-chip \
  -d chip \
  -t cyp3a5 \
  --vcf rok-cyp3a5.vcf
