# variant_file_ETL

This is a python script that extracts data from a vcf file in S3, transforms it to fit a data model, and loads it into an ElasticSearch database.

#### Requirements
pyensembl (EnsemblRelease(75) for hg19, EnsemblRelease(87) for hg38, and EnsemblRelease(54) for hg18), pysam (VariantFile), elasticsearchm and mygene.


#### Summary
The script uses pysams VariantFile to read files from s3. It assumes an index file exists, appended with .tbi. Edit the variable file_list to include all the vcf files for extraction. The script uses Amazons sqs to create a queue of the files. 

Each file is processed one at a time, but rows are processed in parallel. The script determines which genome was used to create the VCF files, parses the info and format fields from the header. It determines the cytoband based from variant position, the gene name from  position, and the protein name from the gene name. This is all in respects to the determined genome used. 

Gene data is pulled from pyensembl (http://useast.ensembl.org/index.html) and protein data is pulled from mygene (https://mygene.info/). 

Finally, it does a bulk upload of all the extracted and transformed data into ElasticScearch for future querying.

