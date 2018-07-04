# An R Shiny server for Bioinformatics apps
A web service accessible through a web browser (from anywhere), that holds a collection of [R shiny](https://shiny.rstudio.com/) apps to perform routine bioinformatics analyses.

## Shinotate
A shiny app to deliver processed and annotated transcriptome data stored in [Trinotate](https://trinotate.github.io/) databases

### Main features
#### Annotation
* Select from available transcriptome databases or custom tables  
* Allow filtering based on keyword/transcript/gene  
* Download selected sequences or annotation table  
* Provide basic summarised information (sequences, annotation, expression)  
#### Expression
* Show selected transcript expression (in a table and graph)  
* Allow consolidating multiple transcripts under one gene name/symbol   
* Provide basic plotting (expression)  
#### Metadata
* Provide basic information on the experimental design (samples, replicates), sequencing platform (read files), etc.  
* Analysis methods/pipeline (Rmarkdown)

See additional information in the [wiki](https://github.com/IdoBar/shiny-server/wiki) page

## CFX-qpcR
A shiny app to analyse qPCR results from Biorad CFX platform.  

### Main features
* Select genes/targets  
* Warn about outliers/extreme variability among replicates  
* Display calibration curves plots  
* Save linear regression, R<sup>2</sup> and calculated *Eff* in a table available for download  
