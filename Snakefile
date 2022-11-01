sample = ["hg002.gangstr"]

rule all:
    input:
        "result/done.txt"

rule sort:
    input :
        "{sample}.bed"
    output :
        temp("{sample}.sorted.bed")
    shell :
        "sort-bed {input} > {output}"


rule anno_gene:
 input :
      "{sample}.sorted.bed"
 output :
      temp("temp/{sample}.1.bed")
 shell:
      "bedmap --bp-ovr 1  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val 'intergenic'  {input} annot_beds/hg38/hg38.UCSCgene.sorted.bed  > {output}"

rule anno_genename:
 input :
      "temp/{sample}.1.bed"
 output :
      temp("temp/{sample}.2.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val '-'  {input} annot_beds/hg38/hg38.UCSCgeneName.sorted.bed  > {output}"

rule anno_CDS:
 input :
      "temp/{sample}.2.bed"
 output :
      temp("temp/{sample}.3.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val '-'  {input} annot_beds/hg38/CDS.hg38.sorted.bed  > {output}"

rule anno_exon:
 input :
      "temp/{sample}.3.bed"
 output :
      temp("temp/{sample}.4.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val '-'  {input} annot_beds/hg38/exon.hg38.sorted.bed  > {output}"

rule anno_UTR:
 input :
      "temp/{sample}.4.bed"
 output :
      temp("temp/{sample}.5.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val '-'  {input} annot_beds/hg38/UTR.hg38.sorted.bed > {output}"

rule anno_ref_motif:
 input :
      "temp/{sample}.5.bed"
 output :
      temp("temp/{sample}.6.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val '-'  {input} annot_beds/hg38/hg38_ver13_GangSTR_STR.motif.sorted.bed > {output}"
      
rule anno_ref_period:
 input :
      "temp/{sample}.6.bed"
 output :
      temp("temp/{sample}.7.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val '-'  {input} annot_beds/hg38/hg38_ver13_GangSTR_STR.period.sorted.bed > {output}"

rule anno_1000g_50Han_mean:
 input :
      "temp/{sample}.7.bed"
 output :
      temp("temp/{sample}.8.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val 'NA'  {input} annot_beds/hg38/1000g_50Han_mean.sorted.bed > {output}"

rule anno_1000g_50Han_sd:
 input :
      "temp/{sample}.8.bed"
 output :
      temp("temp/{sample}.9.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val 'NA'  {input} annot_beds/hg38/1000g_50Han_sd.sorted.bed > {output}"


rule anno_1000g_50Han_max:
 input :
      "temp/{sample}.9.bed"
 output :
      temp("temp/{sample}.10.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val 'NA'  {input} annot_beds/hg38/1000g_50Han_max.sorted.bed > {output}"

rule anno_1000g_100_mean:
 input :
      "temp/{sample}.10.bed"
 output :
      temp("temp/{sample}.11.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val 'NA'  {input} annot_beds/hg38/1000g_100_mean.sorted.bed > {output}"

rule anno_1000g_100_sd:
 input :
      "temp/{sample}.11.bed"
 output :
      temp("temp/{sample}.12.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val 'NA'  {input} annot_beds/hg38/1000g_100_sd.sorted.bed > {output}"


rule anno_1000g_100_max:
 input :
      "temp/{sample}.12.bed"
 output :
      temp("temp/{sample}.13.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val 'NA'  {input} annot_beds/hg38/1000g_100_max.sorted.bed > {output}"


rule anno_TADboundry:
 input :
      "temp/{sample}.13.bed"
 output :
      temp("temp/{sample}.14.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val '-'  {input} annot_beds/hg38/TADBoundaries40kb.sorted.bed > {output}"

rule anno_CTCF:
 input :
      "temp/{sample}.14.bed"
 output :
      temp("temp/{sample}.15.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val '-'  {input} annot_beds/hg38/GM12864CTCF.hg38.sorted.bed > {output}"

rule anno_TFBS:
 input :
      "temp/{sample}.15.bed"
 output :
      temp("temp/{sample}.16.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val '-'  {input} annot_beds/hg38/TFBS.hg38.sorted.bed > {output}"

rule anno_gene_pLI:
 input :
      "temp/{sample}.16.bed"
 output :
      temp("temp/{sample}.17.bed")
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val '-'  {input} annot_beds/hg38/pLI.hg38.sorted.bed > {output}"

rule anno_phenotype:
 input :
      "temp/{sample}.17.bed"
 output :
      "result/{sample}.annot.tsv"
 shell:
      "bedmap  --echo --echo-map-id-uniq  --delim '\t' --unmapped-val '-'  {input} annot_beds/hg38/HPO.hg38.bed > {output}"

rule predict:
 input :
      "result/{sample}.annot.tsv"
 output :
      "result/{sample}.annot.scored.tsv"
 shell:
      "Rscript predict.R {input} {output} "


rule echo_done:
 input :
      expand("result/{sample}.annot.scored.tsv",sample=sample)
 output :
      "result/done.txt"
 shell:
      "echo annoted beds : {input} > {output} "
