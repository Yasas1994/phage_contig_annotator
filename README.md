```
█▀▀█ █░░█ █▀▀█ █▀▀▀ █▀▀ 　 █▀▀ █▀▀█ █▀▀▄ ▀▀█▀▀ ░▀░ █▀▀▀ 　 █▀▀█ █▀▀▄ █▀▀▄ █▀▀█ ▀▀█▀▀ █▀▀█ ▀▀█▀▀ █▀▀█ █▀▀█ 
█░░█ █▀▀█ █▄▄█ █░▀█ █▀▀ 　 █░░ █░░█ █░░█ ░░█░░ ▀█▀ █░▀█ 　 █▄▄█ █░░█ █░░█ █░░█ ░░█░░ █▄▄█ ░░█░░ █░░█ █▄▄▀ 
█▀▀▀ ▀░░▀ ▀░░▀ ▀▀▀▀ ▀▀▀ 　 ▀▀▀ ▀▀▀▀ ▀░░▀ ░░▀░░ ▀▀▀ ▀▀▀▀ 　 ▀░░▀ ▀░░▀ ▀░░▀ ▀▀▀▀ ░░▀░░ ▀░░▀ ░░▀░░ ▀▀▀▀ ▀░▀▀
```
Annotates genes on putative phage contigs with a database of hidden markov models (HMMs) based on [PHROGs](https://phrogs.lmge.uca.fr/). This tool was built to support visual confirmation of predictions made by [Jeager](https://github.com/Yasas1994/Jaeger).

## installation 

1) clone the repository ```git clone https://github.com/Yasas1994/phage_contig_annotator.git && cd phage_contig_annotator```
2) create a conda environment ```conda env create -f environment.yml```
3) install phage_contig_annotator ```pip install .```
4) download database ```phage_contig_annotator download_db```

## usage
```
conda activate phage_contig_annotator 
phage_contig_annotator --cpus 10 runall --contigs --input input_fasta_file_with_contigs --output output_dir

```
## examples 
![image](https://user-images.githubusercontent.com/34155351/201523981-1f5f2fd0-5177-417e-aea0-752430c7eb80.png)
![image](https://user-images.githubusercontent.com/34155351/201524028-13f35b46-83fb-4c50-af81-3701ef992e48.png)
![image](https://user-images.githubusercontent.com/34155351/201524110-c4919fd4-aa34-449b-9ba1-7081a278d6af.png)
