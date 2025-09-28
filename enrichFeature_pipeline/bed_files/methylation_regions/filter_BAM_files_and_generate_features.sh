srcdir="/mnt/DATASM14/DATA_HIEUNGUYEN/src/ecd_wgs_enriched_features";
outputdir="/mnt/DATASM14/DATA_HIEUNGUYEN/outdir/enriched_features_20250928/methylation_regions";
beddir=${srcdir}/enrichFeature_pipeline/bed_files/methylation_regions/processed_methylation_regions
feature_srcdir="/mnt/DATASM14/DATA_HIEUNGUYEN/src/ecd_wgs_feature_pipeline";
path_to_fa="/mnt/archiving/DATA_HIEUNGUYEN/2024/resources/hg19.fa";
path_to_nucleosome_ref="/mnt/archiving/DATA_HIEUNGUYEN/2024/resources/rpr_map_EXP0779.sorted.bed";

metadata=${srcdir}/full_metadata.csv

all_beds=$(ls $beddir/*.bed)

mkdir -p ${outputdir};

files=$(cat ${metadata} | tail -n +2 | cut -d, -f4);

for inputbed in $all_beds;do \
    echo -e "working on bed file " $inputbed;
    bedname=$(echo $inputbed | xargs -n 1 basename);
    mkdir -p ${outputdir}/${bedname};

    for file in $files;do \
        filename=$(echo $file | xargs -n 1 basename);
        echo -e "working on file "  $filename;
        bedname=$(echo $inputbed | xargs -n 1 basename);
        samtools view -b -h -L ${inputbed} ${file} > ${outputdir}/${bedname}/${filename%.bam*}.filtered.bam;
        samtools index ${outputdir}/${bedname}/${filename%.bam*}.filtered.bam;

        ##### generate features
        filtered_bam=${outputdir}/${bedname}/${filename%.bam*}.filtered.bam;
        sampleid=$(echo ${filtered_bam} | xargs -n 1 basename)
        sampleid=$(echo ${sampleid} | cut -d '.' -f 1)

        mkdir -p ${outputdir}/${bedname}/${sampleid};

        echo -e "WORKING ON BAM FILE " $sampleid " with bed file " $bedname;
        bash ${feature_srcdir}/01_convert_BAM_to_FRAGS.sh -i ${filtered_bam} -o ${outputdir}/${bedname} -n 20 -f ${path_to_fa} -q 50 -w 150 -e 151 -r 350 -c "true"

        bash ${feature_srcdir}/02_calculate_features_from_frag.sh -i ${outputdir}/${bedname}/${sampleid}/${sampleid}.frag.tsv -o ${outputdir}/${bedname} -f ${path_to_fa} -r ${path_to_nucleosome_ref} -c "true"

        echo -e "Generating final feature files ..."
        # source /home/hieunguyen/miniconda3/bin/activate && conda activate pytorch;
        echo ${outputdir}/${bedname}/${sampleid}/${sampleid}.frag_output.tsv
        echo  ${outputdir}/${bedname}/${sampleid}/fragmentomics_features
        python ${feature_srcdir}/04_generate_features_fragmentomics_and_image_features.py \
            --input ${outputdir}/${bedname}/${sampleid}/${sampleid}.frag_output.tsv \
            --output ${outputdir}/${bedname}/${sampleid}/fragmentomics_features;
            
        done;done